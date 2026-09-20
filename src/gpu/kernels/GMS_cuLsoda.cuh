/* Copyright (C) Bernard Gingold, 2020-2026 
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#ifndef __GMS_CULSODA_CUH__
#define __GMS_CULSODA_CUH__
#include <math.h>
#include <stdio.h>
#include <cuda_runtime.h>
#include "GMS_gpu_config.cuh"




/*
 
 barf  [ba:rf]  2.  "He suggested using FORTRAN, and everybody barfed."
 
 - From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition) 
 
 */


/* dlsoda.f -- translated by f2c (version 20090411).
 
 http://www.netlib.org/f2c/libf2c.zip
 */


#define double  double


/* Common Block Declarations */
struct cuLsodaCommonBlock
{
	double  /*rowns[209],*/ CM_conit, CM_crate, CM_ccmax, CM_el0, CM_h__, CM_hmin, CM_hmxi, CM_hu, CM_rc, CM_tn, CM_uround, CM_pdest, CM_pdlast, CM_ratio, CM_hold, CM_rmax;
	double   CM_el[13], CM_elco[156]	/* was [13][12] */, CM_tesco[36]	/* was [3][12] */;
	double  CM_rls[218];
	double  CM_tsw, /*rowns2[20],*/ CM_pdnorm;
	double  /*rownd2,*/ CM_cm1[12], CM_cm2[5];
	double  CM_rlsa[22];
	double  CM_sm1[12];
	int CM_init, CM_mxstep, CM_mxhnil, CM_nhnil, CM_nslast, CM_nyh, /*iowns[6],*/ CM_icf, 
	CM_ierpj, CM_iersl, CM_jcur, CM_jstart, CM_kflag, CM_l, CM_lyh, CM_lewt, CM_lacor, CM_lsavf,
	CM_lwm, CM_liwm, CM_meth, CM_miter, CM_maxord, CM_maxcor, CM_msbp, CM_mxncf, CM_n, CM_nq, 
	CM_nst, CM_nfe, CM_nje, CM_nqu;
	int /*iownd[6],*/ CM_ialth, CM_ipup, CM_lmax, /*meo,*/ CM_nqnyh, CM_nslp;
	int CM_ils[37];
	int CM_insufr, CM_insufi, CM_ixpr, /*iowns2[2],*/ CM_jtyp, CM_mused, CM_mxordn, CM_mxords; 
	int /*iownd2[3],*/ CM_icount, CM_irflag;
	int CM_ilsa[9];
};

/* End Common Block */ 


#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))


//#define Fex_and_Jex_definition
//struct myFex
//{
//	__device__ void operator()(int *neq, double  *t, double  *y, double  *ydot/*, void *otherData*/)
//	{
//		ydot[0] = (double )1.0E4 * y[1] * y[2] - (double ).04E0 * y[0];
//		ydot[2] = (double )3.0E7 * y[1] * y[1];
//		ydot[1] = (double )-1.0 * (ydot[0] + ydot[2]);
//	}
//};

//struct myJex
//{
//	__device__ void operator()(int *neq, double  *t, double  *y, int ml, int mu, double  *pd, int nrowpd/*, void *otherData*/)
//	{
//		return;
//	}
//};




/* dlsoda.f -- translated by f2c (version 20090411).
 You must link the resulting object file with libf2c:
 on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm
 or, if you install libf2c.a in a standard place, with -lf2c -lm
 -- in that order, at the end of the command line, as in
 cc *.o -lf2c -lm
 Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,
 
 http://www.netlib.org/f2c/libf2c.zip
 */




/*template<typename Fex, typename Jex>
__device__ int dlsoda_(Fex, int *, double  *, double  *, double  *, int *, double  *, double  *, int *, int *, int *, double  *, int *, int *, int *, Jex, int *, struct cuLsodaCommonBlock *);

template<typename Fex, typename Jex> 
__device__ int dstoda_(int *neq, double  *y, double  *yh, int *NOT_nyh, double  *yh1, double  *ewt, double  *savf, double  *acor, double  *wm, int *iwm, Fex f, Jex jac, struct cuLsodaCommonBlock *common);

template<typename Fex, typename Jex> 
__device__ int dprja_(int *neq, double  *y, double  *yh, int *NOT_nyh, double  *ewt, double  *ftem, double  *savf, double  *wm, int *iwm, Fex f, Jex jac, struct cuLsodaCommonBlock *common);

__device__ int dsolsy_(double  *wm, int *iwm, double  *x, double  *tem, struct cuLsodaCommonBlock *common);
__device__ int dintdy_(double  *t, int k, double  *yh, int *NOT_nyh, double  *dky, int *iflag, struct cuLsodaCommonBlock *common);
__device__ int dcfode_(int meth, double  *DCFODE_elco, double  *DCFODE_tesco, struct cuLsodaCommonBlock *common);
__device__ int dsolsy_(double  *wm, int *iwm, double  *x, double  *tem, struct cuLsodaCommonBlock *common);
__device__ int dewset_(int *PARAM_n, int *itol, double  *rtol, double  *atol, double  *ycur, double  *ewt, struct cuLsodaCommonBlock *common);
__device__ double  dmnorm_(int *PARAM_n, double  *v, double  *w, struct cuLsodaCommonBlock *common);
__device__ double  dfnorm_(int *PARAM_n, double  *a, double  *w, struct cuLsodaCommonBlock *common);
__device__ double  dbnorm_(int *PARAM_n, double  *a, int *nra, int *ml, int *mu, double  *w, struct cuLsodaCommonBlock *common);
__device__ int dsrcma_(double  *rsav, int *isav, int *job, struct cuLsodaCommonBlock *common);
__device__ int dgefa_(double  *a, int *lda, int *PARAM_n, int *ipvt, int *info, struct cuLsodaCommonBlock *common);
__device__ int dgesl_(double  *a, int *lda, int *PARAM_n, int *ipvt, double  *b, int job, struct cuLsodaCommonBlock *common);
__device__ int dgbfa_(double  *abd, int *lda, int *PARAM_n, int *ml, int *mu, int *ipvt, int *info, struct cuLsodaCommonBlock *common);
__device__ int dgbsl_(double  *abd, int *lda, int *PARAM_n, int *ml, int *mu, int *ipvt, double  *b, int job, struct cuLsodaCommonBlock *common);
__device__ double  dumach_( struct cuLsodaCommonBlock *common);
//__device__ int xsetf_(int *mflag, struct CommonBlock *common);
//__device__ int xsetun_(int *lun, struct CommonBlock *common);
__device__ int ixsav_(int ipar, int *ivalue, int iset, struct cuLsodaCommonBlock *common);
__device__ int idamax_(int *PARAM_n, double  *dx, int incx, struct cuLsodaCommonBlock *common);
__device__ int daxpy_(int *PARAM_n, double  *da, double  *dx, int incx, double  *dy, int incy, struct cuLsodaCommonBlock *common);
__device__ int dumsum_(double  a, double  b, double  *c__, struct cuLsodaCommonBlock *common);
__device__ int dscal_(int *PARAM_n, double  *da, double  *dx, int incx, struct cuLsodaCommonBlock *common);
__device__ double  ddot_(int *PARAM_n, double  *dx, int incx, double  *dy, int incy, struct cuLsodaCommonBlock *common);
__device__ double  d_sign(double  *a, double  *b);
__host__ __device__ void cuLsodaCommonBlockInit(struct cuLsodaCommonBlock *common);

*/


 
 
 // Begin Define block for common variables
 #define	conit common->CM_conit
 #define	crate common->CM_crate
 #define	ccmax common->CM_ccmax
 #define	el0 common->CM_el0
 #define	h__ common->CM_h__
 #define	hmin common->CM_hmin
 #define	hmxi common->CM_hmxi
 #define	hu common->CM_hu
 #define	rc common->CM_rc
 #define	tn common->CM_tn
 #define	uround common->CM_uround
 #define	pdest common->CM_pdest
 #define	pdlast common->CM_pdlast
 #define	ratio common->CM_ratio
 #define	hold common->CM_hold
 #define	rmax common->CM_rmax
 #define	el common->CM_el
 #define	elco common->CM_elco
 #define	tesco common->CM_tesco
 #define	rls common->CM_rls
 #define	tsw common->CM_tsw
 #define	pdnorm common->CM_pdnorm
 #define	cm1 common->CM_cm1
 #define	cm2 common->CM_cm2
 #define	rlsa common->CM_rlsa
 #define	sm1 common->CM_sm1
 #define	init common->CM_init
 #define	mxstep common->CM_mxstep
 #define	mxhnil common->CM_mxhnil
 #define	nhnil common->CM_nhnil
 #define	nslast common->CM_nslast
 #define	nyh common->CM_nyh
 #define	icf common->CM_icf
 #define	ierpj common->CM_ierpj
 #define	iersl common->CM_iersl
 #define	jcur common->CM_jcur
 #define	jstart common->CM_jstart
 #define	kflag common->CM_kflag
 #define	l common->CM_l
 #define	lyh common->CM_lyh
 #define	lewt common->CM_lewt
 #define	lacor common->CM_lacor
 #define	lsavf common->CM_lsavf
 #define	lwm common->CM_lwm
 #define	liwm common->CM_liwm
 #define	meth common->CM_meth
 #define	miter common->CM_miter
 #define	maxord common->CM_maxord
 #define	maxcor common->CM_maxcor
 #define	msbp common->CM_msbp
 #define	mxncf common->CM_mxncf
 #define	n common->CM_n
 #define	nq common->CM_nq
 #define	nst common->CM_nst
 #define	nfe common->CM_nfe
 #define	nje common->CM_nje
 #define	nqu common->CM_nqu
 #define	ialth common->CM_ialth
 #define	ipup common->CM_ipup
 #define	lmax common->CM_lmax
 #define	nqnyh common->CM_nqnyh
 #define	nslp common->CM_nslp
 #define	ils common->CM_ils
 #define	insufr common->CM_insufr
 #define	insufi common->CM_insufi
 #define	ixpr common->CM_ixpr
 #define	jtyp common->CM_jtyp
 #define	mused common->CM_mused
 #define	mxordn common->CM_mxordn
 #define	mxords common->CM_mxords
 #define	icount common->CM_icount
 #define	irflag common->CM_irflag
 #define	ilsa common->CM_ilsa
 // End of Definitions
 


    //struct myFex
//{
//	__device__ void operator()(int *neq, double  *t, double  *y, double  *ydot/*, void *otherData*/)
//	{
//		ydot[0] = (double )1.0E4 * y[1] * y[2] - (double ).04E0 * y[0];
//		ydot[2] = (double )3.0E7 * y[1] * y[1];
//		ydot[1] = (double )-1.0 * (ydot[0] + ydot[2]);
//	}
//};

//struct myJex
//{
//	__device__ void operator()(int *neq, double  *t, double  *y, int ml, int mu, double  *pd, int nrowpd/*, void *otherData*/)
//	{
//		return;
//	}
//};




template<typename Fex,typename Jex>
__device__ int dlsoda_(Fex f, 
                       int *  __restrict__ neq, 
                       double  * __restrict__ y, 
                       double  * __restrict__ t, 
                       double  * __restrict__ tout, 
                       int *  __restrict__ itol, 
                       double  * __restrict__ rtol, 
                       double  * __restrict__ atol, 
                       int * __restrict__ itask, 
                       int * __restrict__ istate, 
                       int * __restrict__ iopt, 
                       double  * __restrict__ rwork, 
                       int * __restrict__ lrw, 
                       int * __restrict__ iwork, 
                       int * __restrict__ liw, 
                       Jex jac, 
                       int * __restrict__ jt, 
                       struct cuLsodaCommonBlock * __restrict__ common)
{
    /* Initialized data */
	//struct cuLsodaCommonBlock commonB;
	//struct cuLsodaCommonBlock *common;
	//common = &commonB;
		
	int errorCode = 0;
	
	int mord[2] = { 12,5 };
	int mxstp0 = 500;
	int mxhnl0 = 10;
	
    /* System generated locals */
    int i__1 = 0;
    double  d__1 = 0.;
	double  d__2 = 0.;

	
	
    /* Local variables */
     int i__;
     double  h0 = 0.;
     int i1 = 0;
	 int i2 = 0;
     double  w0 = 0.;
     int ml = 0;
     double  rh = 0.;
     int mu = 0;
     double  tp = 0.;
     int lf0 = 0;
     double  big = 0.;
     int kgo = 0;
     double  ayi = 0.;
     double  hmx = 0.;
	 double  tol = 0.;
	 double  sum = 0.;
     int len1 = 0;
	 int len2 = 0;
     double  hmax = 0.;
     int ihit = 0;
     double  ewti = 0.;
	 double  size = 0.;
     int len1c = 0;
	 int len1n = 0;
	 int len1s = 0;
	 int iflag;
     double  atoli = 0.;
     int leniw = 0;
	 int lenwm = 0;
	 int imxer = 0;
     double  tcrit = 0.;
     int lenrw = 0;
     double  tdist = 0.;
	 double  rtoli = 0.; 
	 double  tolsf = 0.;
	 double  tnext = 0.;
	
     int leniwc = 0;
     int lenrwc = 0;
	
	/* ----------------------------------------------------------------------- */
	/* This is the 12 November 2003 version of */
	/* DLSODA: Livermore Solver for Ordinary Differential Equations, with */
	/*         Automatic method switching for stiff and nonstiff problems. */
	
	/* This version is in double  precision. */
	
	/* DLSODA solves the initial value problem for stiff or nonstiff */
	/* systems of first order ODEs, */
	/*     dy/dt = f(t,y) ,  or, in component form, */
	/*     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ). */
	
	/* This a variant version of the DLSODE package. */
	/* It switches automatically between stiff and nonstiff methods. */
	/* This means that the user does not have to determine whether the */
	/* problem is stiff or not, and the solver will automatically choose the */
	/* appropriate method.  It always starts with the nonstiff method. */
	
	/* Authors:       Alan C. Hindmarsh */
	/*                Center for Applied Scientific Computing, L-561 */
	/*                Lawrence Livermore National Laboratory */
	/*                Livermore, CA 94551 */
	/* and */
	/*                Linda R. Petzold */
	/*                Univ. of California at Santa Barbara */
	/*                Dept. of Computer Science */
	/*                Santa Barbara, CA 93106 */
	
	/* References: */
	/* 1.  Alan C. Hindmarsh,  ODEPACK, A Systematized Collection of ODE */
	/*     Solvers, in Scientific Computing, R. S. Stepleman et al. (Eds.), */
	/*     North-Holland, Amsterdam, 1983, pp. 55-64. */
	/* 2.  Linda R. Petzold, Automatic Selection of Methods for Solving */
	/*     Stiff and Nonstiff Systems of Ordinary Differential Equations, */
	/*     Siam J. Sci. Stat. Comput. 4 (1983), pp. 136-148. */
	/* ----------------------------------------------------------------------- */
	/* Summary of Usage. */
	
	/* Communication between the user and the DLSODA package, for normal */
	/* situations, is summarized here.  This summary describes only a subset */
	/* of the full set of options available.  See the full description for */
	/* details, including alternative treatment of the Jacobian matrix, */
	/* optional inputs and outputs, nonstandard options, and */
	/* instructions for special situations.  See also the example */
	/* problem (with program and output) following this summary. */
	
	/* A. First provide a subroutine of the form: */
	/*               SUBROUTINE F (NEQ, T, Y, YDOT) */
	/*               DOUBLE PRECISION T, Y(*), YDOT(*) */
	/* which supplies the vector function f by loading YDOT(i) with f(i). */
	
	/* B. Write a main program which calls Subroutine DLSODA once for */
	/* each point at which answers are desired.  This should also provide */
	/* for possible use of logical unit 6 for output of error messages */
	/* by DLSODA.  On the first call to DLSODA, supply arguments as follows: */
	/* F      = name of subroutine for right-hand side vector f. */
	/*          This name must be declared External in calling program. */
	/* NEQ    = number of first order ODEs. */
	/* Y      = array of initial values, of length NEQ. */
	/* T      = the initial value of the independent variable. */
	/* TOUT   = first point where output is desired (.ne. T). */
	/* ITOL   = 1 or 2 according as ATOL (below) is a scalar or array. */
	/* RTOL   = relative tolerance parameter (scalar). */
	/* ATOL   = absolute tolerance parameter (scalar or array). */
	/*          the estimated local error in y(i) will be controlled so as */
	/*          to be less than */
	/*             EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or */
	/*             EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2. */
	/*          Thus the local error test passes if, in each component, */
	/*          either the absolute error is less than ATOL (or ATOL(i)), */
	/*          or the relative error is less than RTOL. */
	/*          Use RTOL = 0.0 for pure absolute error control, and */
	/*          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error */
	/*          control.  Caution: actual (global) errors may exceed these */
	/*          local tolerances, so choose them conservatively. */
	/* ITASK  = 1 for normal computation of output values of y at t = TOUT. */
	/* ISTATE = int flag (input and output).  Set ISTATE = 1. */
	/* IOPT   = 0 to indicate no optional inputs used. */
	/* RWORK  = real work array of length at least: */
	/*             22 + NEQ * MAX(16, NEQ + 9). */
	/*          See also Paragraph E below. */
	/* LRW    = declared length of RWORK (in user's dimension). */
	/* IWORK  = int work array of length at least  20 + NEQ. */
	/* LIW    = declared length of IWORK (in user's dimension). */
	/* JAC    = name of subroutine for Jacobian matrix. */
	/*          Use a dummy name.  See also Paragraph E below. */
	/* JT     = Jacobian type indicator.  Set JT = 2. */
	/*          See also Paragraph E below. */
	/* Note that the main program must declare arrays Y, RWORK, IWORK, */
	/* and possibly ATOL. */
	
	/* C. The output from the first call (or any call) is: */
	/*      Y = array of computed values of y(t) vector. */
	/*      T = corresponding value of independent variable (normally TOUT). */
	/* ISTATE = 2  if DLSODA was successful, negative otherwise. */
	/*          -1 means excess work done on this call (perhaps wrong JT). */
	/*          -2 means excess accuracy requested (tolerances too small). */
	/*          -3 means illegal input detected (see printed message). */
	/*          -4 means repeated error test failures (check all inputs). */
	/*          -5 means repeated convergence failures (perhaps bad Jacobian */
	/*             supplied or wrong choice of JT or tolerances). */
	/*          -6 means error weight became zero during problem. (Solution */
	/*             component i vanished, and ATOL or ATOL(i) = 0.) */
	/*          -7 means work space insufficient to finish (see messages). */
	
	/* D. To continue the integration after a successful return, simply */
	/* reset TOUT and call DLSODA again.  No other parameters need be reset. */
	
	/* E. Note: If and when DLSODA regards the problem as stiff, and */
	/* switches methods accordingly, it must make use of the NEQ by NEQ */
	/* Jacobian matrix, J = df/dy.  For the sake of simplicity, the */
	/* inputs to DLSODA recommended in Paragraph B above cause DLSODA to */
	/* treat J as a full matrix, and to approximate it internally by */
	/* difference quotients.  Alternatively, J can be treated as a band */
	/* matrix (with great potential reduction in the size of the RWORK */
	/* array).  Also, in either the full or banded case, the user can supply */
	/* J in closed form, with a routine whose name is passed as the JAC */
	/* argument.  These alternatives are described in the paragraphs on */
	/* RWORK, JAC, and JT in the full description of the call sequence below. */
	
	/* ----------------------------------------------------------------------- */
	/* Example Problem. */
	
	/* The following is a simple example problem, with the coding */
	/* needed for its solution by DLSODA.  The problem is from chemical */
	/* kinetics, and consists of the following three rate equations: */
	/*     dy1/dt = -.04*y1 + 1.e4*y2*y3 */
	/*     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2 */
	/*     dy3/dt = 3.e7*y2**2 */
	/* on the interval from t = 0.0 to t = 4.e10, with initial conditions */
	/* y1 = 1.0, y2 = y3 = 0.  The problem is stiff. */
	
	/* The following coding solves this problem with DLSODA, */
	/* printing results at t = .4, 4., ..., 4.e10.  It uses */
	/* ITOL = 2 and ATOL much smaller for y2 than y1 or y3 because */
	/* y2 has much smaller values. */
	/* At the end of the run, statistical quantities of interest are */
	/* printed (see optional outputs in the full description below). */
	
	/*     EXTERNAL FEX */
	/*     DOUBLE PRECISION ATOL, RTOL, RWORK, T, TOUT, Y */
	/*     DIMENSION Y(3), ATOL(3), RWORK(70), IWORK(23) */
	/*     NEQ = 3 */
	/*     Y(1) = 1. */
	/*     Y(2) = 0. */
	/*     Y(3) = 0. */
	/*     T = 0. */
	/*     TOUT = .4 */
	/*     ITOL = 2 */
	/*     RTOL = 1.D-4 */
	/*     ATOL(1) = 1.D-6 */
	/*     ATOL(2) = 1.D-10 */
	/*     ATOL(3) = 1.D-6 */
	/*     ITASK = 1 */
	/*     ISTATE = 1 */
	/*     IOPT = 0 */
	/*     LRW = 70 */
	/*     LIW = 23 */
	/*     JT = 2 */
	/*     DO 40 IOUT = 1,12 */
	/*       CALL DLSODA(FEX,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE, */
	/*    1     IOPT,RWORK,LRW,IWORK,LIW,JDUM,JT) */
	/*       WRITE(6,20)T,Y(1),Y(2),Y(3) */
	/* 20    FORMAT(' At t =',D12.4,'   Y =',3D14.6) */
	/*       IF (ISTATE .LT. 0) GO TO 80 */
	/* 40    TOUT = TOUT*10. */
	/*     WRITE(6,60)IWORK(11),IWORK(12),IWORK(13),IWORK(19),RWORK(15) */
	/* 60  FORMAT(/' No. steps =',I4,'  No. f-s =',I4,'  No. J-s =',I4/ */
	/*    1   ' Method last used =',I2,'   Last switch was at t =',D12.4) */
	/*     STOP */
	/* 80  WRITE(6,90)ISTATE */
	/* 90  FORMAT(///' Error halt.. ISTATE =',I3) */
	/*     STOP */
	/*     END */
	
	/*     SUBROUTINE FEX (NEQ, T, Y, YDOT) */
	/*     DOUBLE PRECISION T, Y, YDOT */
	/*     DIMENSION Y(3), YDOT(3) */
	/*     YDOT(1) = -.04*Y(1) + 1.D4*Y(2)*Y(3) */
	/*     YDOT(3) = 3.D7*Y(2)*Y(2) */
	/*     YDOT(2) = -YDOT(1) - YDOT(3) */
	/*     RETURN */
	/*     END */
	
	/* The output of this program (on a CDC-7600 in single precision) */
	/* is as follows: */
	
	/*   At t =  4.0000e-01   y =  9.851712e-01  3.386380e-05  1.479493e-02 */
	/*   At t =  4.0000e+00   Y =  9.055333e-01  2.240655e-05  9.444430e-02 */
	/*   At t =  4.0000e+01   Y =  7.158403e-01  9.186334e-06  2.841505e-01 */
	/*   At t =  4.0000e+02   Y =  4.505250e-01  3.222964e-06  5.494717e-01 */
	/*   At t =  4.0000e+03   Y =  1.831975e-01  8.941774e-07  8.168016e-01 */
	/*   At t =  4.0000e+04   Y =  3.898730e-02  1.621940e-07  9.610125e-01 */
	/*   At t =  4.0000e+05   Y =  4.936363e-03  1.984221e-08  9.950636e-01 */
	/*   At t =  4.0000e+06   Y =  5.161831e-04  2.065786e-09  9.994838e-01 */
	/*   At t =  4.0000e+07   Y =  5.179817e-05  2.072032e-10  9.999482e-01 */
	/*   At t =  4.0000e+08   Y =  5.283401e-06  2.113371e-11  9.999947e-01 */
	/*   At t =  4.0000e+09   Y =  4.659031e-07  1.863613e-12  9.999995e-01 */
	/*   At t =  4.0000e+10   Y =  1.404280e-08  5.617126e-14  1.000000e+00 */
	
	/*   No. steps = 361  No. f-s = 693  No. J-s =  64 */
	/*   Method last used = 2   Last switch was at t =  6.0092e-03 */
	/* ----------------------------------------------------------------------- */
	/* Full description of user interface to DLSODA. */
	
	/* The user interface to DLSODA consists of the following parts. */
	
	/* 1.   The call sequence to Subroutine DLSODA, which is a driver */
	/*      routine for the solver.  This includes descriptions of both */
	/*      the call sequence arguments and of user-supplied routines. */
	/*      following these descriptions is a description of */
	/*      optional inputs available through the call sequence, and then */
	/*      a description of optional outputs (in the work arrays). */
	
	/* 2.   Descriptions of other routines in the DLSODA package that may be */
	/*      (optionally) called by the user.  These provide the ability to */
	/*      alter error message handling, save and restore the internal */
	/*      Common, and obtain specified derivatives of the solution y(t). */
	
	/* 3.   Descriptions of Common blocks to be declared in overlay */
	/*      or similar environments, or to be saved when doing an interrupt */
	/*      of the problem and continued solution later. */
	
	/* 4.   Description of a subroutine in the DLSODA package, */
	/*      which the user may replace with his/her own version, if desired. */
	/*      this relates to the measurement of errors. */
	
	/* ----------------------------------------------------------------------- */
	/* Part 1.  Call Sequence. */
	
	/* The call sequence parameters used for input only are */
	/*     F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC, JT, */
	/* and those used for both input and output are */
	/*     Y, T, ISTATE. */
	/* The work arrays RWORK and IWORK are also used for conditional and */
	/* optional inputs and optional outputs.  (The term output here refers */
	/* to the return from Subroutine DLSODA to the user's calling program.) */
	
	/* The legality of input parameters will be thoroughly checked on the */
	/* initial call for the problem, but not checked thereafter unless a */
	/* change in input parameters is flagged by ISTATE = 3 on input. */
	
	/* The descriptions of the call arguments are as follows. */
	
	/* F      = the name of the user-supplied subroutine defining the */
	/*          ODE system.  The system must be put in the first-order */
	/*          form dy/dt = f(t,y), where f is a vector-valued function */
	/*          of the scalar t and the vector y.  Subroutine F is to */
	/*          compute the function f.  It is to have the form */
	/*               SUBROUTINE F (NEQ, T, Y, YDOT) */
	/*               DOUBLE PRECISION T, Y(*), YDOT(*) */
	/*          where NEQ, T, and Y are input, and the array YDOT = f(t,y) */
	/*          is output.  Y and YDOT are arrays of length NEQ. */
	/*          Subroutine F should not alter Y(1),...,Y(NEQ). */
	/*          F must be declared External in the calling program. */
	
	/*          Subroutine F may access user-defined quantities in */
	/*          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array */
	/*          (dimensioned in F) and/or Y has length exceeding NEQ(1). */
	/*          See the descriptions of NEQ and Y below. */
	
	/*          If quantities computed in the F routine are needed */
	/*          externally to DLSODA, an extra call to F should be made */
	/*          for this purpose, for consistent and accurate results. */
	/*          If only the derivative dy/dt is needed, use DINTDY instead. */
	
	/* NEQ    = the size of the ODE system (number of first order */
	/*          ordinary differential equations).  Used only for input. */
	/*          NEQ may be decreased, but not increased, during the problem. */
	/*          If NEQ is decreased (with ISTATE = 3 on input), the */
	/*          remaining components of Y should be left undisturbed, if */
	/*          these are to be accessed in F and/or JAC. */
	
	/*          Normally, NEQ is a scalar, and it is generally referred to */
	/*          as a scalar in this user interface description.  However, */
	/*          NEQ may be an array, with NEQ(1) set to the system size. */
	/*          (The DLSODA package accesses only NEQ(1).)  In either case, */
	/*          this parameter is passed as the NEQ argument in all calls */
	/*          to F and JAC.  Hence, if it is an array, locations */
	/*          NEQ(2),... may be used to store other int data and pass */
	/*          it to F and/or JAC.  Subroutines F and/or JAC must include */
	/*          NEQ in a Dimension statement in that case. */
	
	/* Y      = a real array for the vector of dependent variables, of */
	/*          length NEQ or more.  Used for both input and output on the */
	/*          first call (ISTATE = 1), and only for output on other calls. */
	/*          On the first call, Y must contain the vector of initial */
	/*          values.  On output, Y contains the computed solution vector, */
	/*          evaluated at T.  If desired, the Y array may be used */
	/*          for other purposes between calls to the solver. */
	
	/*          This array is passed as the Y argument in all calls to */
	/*          F and JAC.  Hence its length may exceed NEQ, and locations */
	/*          Y(NEQ+1),... may be used to store other real data and */
	/*          pass it to F and/or JAC.  (The DLSODA package accesses only */
	/*          Y(1),...,Y(NEQ).) */
	
	/* T      = the independent variable.  On input, T is used only on the */
	/*          first call, as the initial point of the integration. */
	/*          on output, after each call, T is the value at which a */
	/*          computed solution Y is evaluated (usually the same as TOUT). */
	/*          on an error return, T is the farthest point reached. */
	
	/* TOUT   = the next value of t at which a computed solution is desired. */
	/*          Used only for input. */
	
	/*          When starting the problem (ISTATE = 1), TOUT may be equal */
	/*          to T for one call, then should .ne. T for the next call. */
	/*          For the initial t, an input value of TOUT .ne. T is used */
	/*          in order to determine the direction of the integration */
	/*          (i.e. the algebraic sign of the step sizes) and the rough */
	/*          scale of the problem.  Integration in either direction */
	/*          (forward or backward in t) is permitted. */
	
	/*          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after */
	/*          the first call (i.e. the first call with TOUT .ne. T). */
	/*          Otherwise, TOUT is required on every call. */
	
	/*          If ITASK = 1, 3, or 4, the values of TOUT need not be */
	/*          monotone, but a value of TOUT which backs up is limited */
	/*          to the current internal T interval, whose endpoints are */
	/*          TCUR - HU and TCUR (see optional outputs, below, for */
	/*          TCUR and HU). */
	
	/* ITOL   = an indicator for the type of error control.  See */
	/*          description below under ATOL.  Used only for input. */
	
	/* RTOL   = a relative error tolerance parameter, either a scalar or */
	/*          an array of length NEQ.  See description below under ATOL. */
	/*          Input only. */
	
	/* ATOL   = an absolute error tolerance parameter, either a scalar or */
	/*          an array of length NEQ.  Input only. */
	
	/*             The input parameters ITOL, RTOL, and ATOL determine */
	/*          the error control performed by the solver.  The solver will */
	/*          control the vector E = (E(i)) of estimated local errors */
	/*          in y, according to an inequality of the form */
	/*                      max-norm of ( E(i)/EWT(i) )   .le.   1, */
	/*          where EWT = (EWT(i)) is a vector of positive error weights. */
	/*          The values of RTOL and ATOL should all be non-negative. */
	/*          The following table gives the types (scalar/array) of */
	/*          RTOL and ATOL, and the corresponding form of EWT(i). */
	
	/*             ITOL    RTOL       ATOL          EWT(i) */
	/*              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL */
	/*              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i) */
	/*              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL */
	/*              4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i) */
	
	/*          When either of these parameters is a scalar, it need not */
	/*          be dimensioned in the user's calling program. */
	
	/*          If none of the above choices (with ITOL, RTOL, and ATOL */
	/*          fixed throughout the problem) is suitable, more general */
	/*          error controls can be obtained by substituting a */
	/*          user-supplied routine for the setting of EWT. */
	/*          See Part 4 below. */
	
	/*          If global errors are to be estimated by making a repeated */
	/*          run on the same problem with smaller tolerances, then all */
	/*          components of RTOL and ATOL (i.e. of EWT) should be scaled */
	/*          down uniformly. */
	
	/* ITASK  = an index specifying the task to be performed. */
	/*          Input only.  ITASK has the following values and meanings. */
	/*          1  means normal computation of output values of y(t) at */
	/*             t = TOUT (by overshooting and interpolating). */
	/*          2  means take one step only and return. */
	/*          3  means stop at the first internal mesh point at or */
	/*             beyond t = TOUT and return. */
	/*          4  means normal computation of output values of y(t) at */
	/*             t = TOUT but without overshooting t = TCRIT. */
	/*             TCRIT must be input as RWORK(1).  TCRIT may be equal to */
	/*             or beyond TOUT, but not behind it in the direction of */
	/*             integration.  This option is useful if the problem */
	/*             has a singularity at or beyond t = TCRIT. */
	/*          5  means take one step, without passing TCRIT, and return. */
	/*             TCRIT must be input as RWORK(1). */
	
	/*          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT */
	/*          (within roundoff), it will return T = TCRIT (exactly) to */
	/*          indicate this (unless ITASK = 4 and TOUT comes before TCRIT, */
	/*          in which case answers at t = TOUT are returned first). */
	
	/* ISTATE = an index used for input and output to specify the */
	/*          the state of the calculation. */
	
	/*          On input, the values of ISTATE are as follows. */
	/*          1  means this is the first call for the problem */
	/*             (initializations will be done).  See note below. */
	/*          2  means this is not the first call, and the calculation */
	/*             is to continue normally, with no change in any input */
	/*             parameters except possibly TOUT and ITASK. */
	/*             (If ITOL, RTOL, and/or ATOL are changed between calls */
	/*             with ISTATE = 2, the new values will be used but not */
	/*             tested for legality.) */
	/*          3  means this is not the first call, and the */
	/*             calculation is to continue normally, but with */
	/*             a change in input parameters other than */
	/*             TOUT and ITASK.  Changes are allowed in */
	/*             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, JT, ML, MU, */
	/*             and any optional inputs except H0, MXORDN, and MXORDS. */
	/*             (See IWORK description for ML and MU.) */
	/*          Note:  A preliminary call with TOUT = T is not counted */
	/*          as a first call here, as no initialization or checking of */
	/*          input is done.  (Such a call is sometimes useful for the */
	/*          purpose of outputting the initial conditions.) */
	/*          Thus the first call for which TOUT .ne. T requires */
	/*          ISTATE = 1 on input. */
	
	/*          On output, ISTATE has the following values and meanings. */
	/*           1  means nothing was done; TOUT = T and ISTATE = 1 on input. */
	/*           2  means the integration was performed successfully. */
	/*          -1  means an excessive amount of work (more than MXSTEP */
	/*              steps) was done on this call, before completing the */
	/*              requested task, but the integration was otherwise */
	/*              successful as far as T.  (MXSTEP is an optional input */
	/*              and is normally 500.)  To continue, the user may */
	/*              simply reset ISTATE to a value .gt. 1 and call again */
	/*              (the excess work step counter will be reset to 0). */
	/*              In addition, the user may increase MXSTEP to avoid */
	/*              this error return (see below on optional inputs). */
	/*          -2  means too much accuracy was requested for the precision */
	/*              of the machine being used.  This was detected before */
	/*              completing the requested task, but the integration */
	/*              was successful as far as T.  To continue, the tolerance */
	/*              parameters must be reset, and ISTATE must be set */
	/*              to 3.  The optional output TOLSF may be used for this */
	/*              purpose.  (Note: If this condition is detected before */
	/*              taking any steps, then an illegal input return */
	/*              (ISTATE = -3) occurs instead.) */
	/*          -3  means illegal input was detected, before taking any */
	/*              integration steps.  See written message for details. */
	/*              Note:  If the solver detects an infinite loop of calls */
	/*              to the solver with illegal input, it will cause */
	/*              the run to stop. */
	/*          -4  means there were repeated error test failures on */
	/*              one attempted step, before completing the requested */
	/*              task, but the integration was successful as far as T. */
	/*              The problem may have a singularity, or the input */
	/*              may be inappropriate. */
	/*          -5  means there were repeated convergence test failures on */
	/*              one attempted step, before completing the requested */
	/*              task, but the integration was successful as far as T. */
	/*              This may be caused by an inaccurate Jacobian matrix, */
	/*              if one is being used. */
	/*          -6  means EWT(i) became zero for some i during the */
	/*              integration.  Pure relative error control (ATOL(i)=0.0) */
	/*              was requested on a variable which has now vanished. */
	/*              The integration was successful as far as T. */
	/*          -7  means the length of RWORK and/or IWORK was too small to */
	/*              proceed, but the integration was successful as far as T. */
	/*              This happens when DLSODA chooses to switch methods */
	/*              but LRW and/or LIW is too small for the new method. */
	
	/*          Note:  Since the normal output value of ISTATE is 2, */
	/*          it does not need to be reset for normal continuation. */
	/*          Also, since a negative input value of ISTATE will be */
	/*          regarded as illegal, a negative output value requires the */
	/*          user to change it, and possibly other inputs, before */
	/*          calling the solver again. */
	
	/* IOPT   = an int flag to specify whether or not any optional */
	/*          inputs are being used on this call.  Input only. */
	/*          The optional inputs are listed separately below. */
	/*          IOPT = 0 means no optional inputs are being used. */
	/*                   default values will be used in all cases. */
	/*          IOPT = 1 means one or more optional inputs are being used. */
	
	/* RWORK  = a real array (double  precision) for work space, and (in the */
	/*          first 20 words) for conditional and optional inputs and */
	/*          optional outputs. */
	/*          As DLSODA switches automatically between stiff and nonstiff */
	/*          methods, the required length of RWORK can change during the */
	/*          problem.  Thus the RWORK array passed to DLSODA can either */
	/*          have a static (fixed) length large enough for both methods, */
	/*          or have a dynamic (changing) length altered by the calling */
	/*          program in response to output from DLSODA. */
	
	/*                       --- Fixed Length Case --- */
	/*          If the RWORK length is to be fixed, it should be at least */
	/*               MAX (LRN, LRS), */
	/*          where LRN and LRS are the RWORK lengths required when the */
	/*          current method is nonstiff or stiff, respectively. */
	
	/*          The separate RWORK length requirements LRN and LRS are */
	/*          as follows: */
	/*          IF NEQ is constant and the maximum method orders have */
	/*          their default values, then */
	/*             LRN = 20 + 16*NEQ, */
	/*             LRS = 22 + 9*NEQ + NEQ**2           if JT = 1 or 2, */
	/*             LRS = 22 + 10*NEQ + (2*ML+MU)*NEQ   if JT = 4 or 5. */
	/*          Under any other conditions, LRN and LRS are given by: */
	/*             LRN = 20 + NYH*(MXORDN+1) + 3*NEQ, */
	/*             LRS = 20 + NYH*(MXORDS+1) + 3*NEQ + LMAT, */
	/*          where */
	/*             NYH    = the initial value of NEQ, */
	/*             MXORDN = 12, unless a smaller value is given as an */
	/*                      optional input, */
	/*             MXORDS = 5, unless a smaller value is given as an */
	/*                      optional input, */
	/*             LMAT   = length of matrix work space: */
	/*             LMAT   = NEQ**2 + 2              if JT = 1 or 2, */
	/*             LMAT   = (2*ML + MU + 1)*NEQ + 2 if JT = 4 or 5. */
	
	/*                       --- Dynamic Length Case --- */
	/*          If the length of RWORK is to be dynamic, then it should */
	/*          be at least LRN or LRS, as defined above, depending on the */
	/*          current method.  Initially, it must be at least LRN (since */
	/*          DLSODA starts with the nonstiff method).  On any return */
	/*          from DLSODA, the optional output MCUR indicates the current */
	/*          method.  If MCUR differs from the value it had on the */
	/*          previous return, or if there has only been one call to */
	/*          DLSODA and MCUR is now 2, then DLSODA has switched */
	/*          methods during the last call, and the length of RWORK */
	/*          should be reset (to LRN if MCUR = 1, or to LRS if */
	/*          MCUR = 2).  (An increase in the RWORK length is required */
	/*          if DLSODA returned ISTATE = -7, but not otherwise.) */
	/*          After resetting the length, call DLSODA with ISTATE = 3 */
	/*          to signal that change. */
	
	/* LRW    = the length of the array RWORK, as declared by the user. */
	/*          (This will be checked by the solver.) */
	
	/* IWORK  = an int array for work space. */
	/*          As DLSODA switches automatically between stiff and nonstiff */
	/*          methods, the required length of IWORK can change during */
	/*          problem, between */
	/*             LIS = 20 + NEQ   and   LIN = 20, */
	/*          respectively.  Thus the IWORK array passed to DLSODA can */
	/*          either have a fixed length of at least 20 + NEQ, or have a */
	/*          dynamic length of at least LIN or LIS, depending on the */
	/*          current method.  The comments on dynamic length under */
	/*          RWORK above apply here.  Initially, this length need */
	/*          only be at least LIN = 20. */
	
	/*          The first few words of IWORK are used for conditional and */
	/*          optional inputs and optional outputs. */
	
	/*          The following 2 words in IWORK are conditional inputs: */
	/*            IWORK(1) = ML     these are the lower and upper */
	/*            IWORK(2) = MU     half-bandwidths, respectively, of the */
	/*                       banded Jacobian, excluding the main diagonal. */
	/*                       The band is defined by the matrix locations */
	/*                       (i,j) with i-ML .le. j .le. i+MU.  ML and MU */
	/*                       must satisfy  0 .le.  ML,MU  .le. NEQ-1. */
	/*                       These are required if JT is 4 or 5, and */
	/*                       ignored otherwise.  ML and MU may in fact be */
	/*                       the band parameters for a matrix to which */
	/*                       df/dy is only approximately equal. */
	
	/* LIW    = the length of the array IWORK, as declared by the user. */
	/*          (This will be checked by the solver.) */
	
	/* Note: The base addresses of the work arrays must not be */
	/* altered between calls to DLSODA for the same problem. */
	/* The contents of the work arrays must not be altered */
	/* between calls, except possibly for the conditional and */
	/* optional inputs, and except for the last 3*NEQ words of RWORK. */
	/* The latter space is used for internal scratch space, and so is */
	/* available for use by the user outside DLSODA between calls, if */
	/* desired (but not for use by F or JAC). */
	
	/* JAC    = the name of the user-supplied routine to compute the */
	/*          Jacobian matrix, df/dy, if JT = 1 or 4.  The JAC routine */
	/*          is optional, but if the problem is expected to be stiff much */
	/*          of the time, you are encouraged to supply JAC, for the sake */
	/*          of efficiency.  (Alternatively, set JT = 2 or 5 to have */
	/*          DLSODA compute df/dy internally by difference quotients.) */
	/*          If and when DLSODA uses df/dy, it treats this NEQ by NEQ */
	/*          matrix either as full (JT = 1 or 2), or as banded (JT = */
	/*          4 or 5) with half-bandwidths ML and MU (discussed under */
	/*          IWORK above).  In either case, if JT = 1 or 4, the JAC */
	/*          routine must compute df/dy as a function of the scalar t */
	/*          and the vector y.  It is to have the form */
	/*               SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD) */
	/*               DOUBLE PRECISION T, Y(*), PD(NROWPD,*) */
	/*          where NEQ, T, Y, ML, MU, and NROWPD are input and the array */
	/*          PD is to be loaded with partial derivatives (elements of */
	/*          the Jacobian matrix) on output.  PD must be given a first */
	/*          dimension of NROWPD.  T and Y have the same meaning as in */
	/*          Subroutine F. */
	/*               In the full matrix case (JT = 1), ML and MU are */
	/*          ignored, and the Jacobian is to be loaded into PD in */
	/*          columnwise manner, with df(i)/dy(j) loaded into PD(i,j). */
	/*               In the band matrix case (JT = 4), the elements */
	/*          within the band are to be loaded into PD in columnwise */
	/*          manner, with diagonal lines of df/dy loaded into the rows */
	/*          of PD.  Thus df(i)/dy(j) is to be loaded into PD(i-j+MU+1,j). */
	/*          ML and MU are the half-bandwidth parameters (see IWORK). */
	/*          The locations in PD in the two triangular areas which */
	/*          correspond to nonexistent matrix elements can be ignored */
	/*          or loaded arbitrarily, as they are overwritten by DLSODA. */
	/*               JAC need not provide df/dy exactly.  A crude */
	/*          approximation (possibly with a smaller bandwidth) will do. */
	/*               In either case, PD is preset to zero by the solver, */
	/*          so that only the nonzero elements need be loaded by JAC. */
	/*          Each call to JAC is preceded by a call to F with the same */
	/*          arguments NEQ, T, and Y.  Thus to gain some efficiency, */
	/*          intermediate quantities shared by both calculations may be */
	/*          saved in a user Common block by F and not recomputed by JAC, */
	/*          if desired.  Also, JAC may alter the Y array, if desired. */
	/*          JAC must be declared External in the calling program. */
	/*               Subroutine JAC may access user-defined quantities in */
	/*          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array */
	/*          (dimensioned in JAC) and/or Y has length exceeding NEQ(1). */
	/*          See the descriptions of NEQ and Y above. */
	
	/* JT     = Jacobian type indicator.  Used only for input. */
	/*          JT specifies how the Jacobian matrix df/dy will be */
	/*          treated, if and when DLSODA requires this matrix. */
	/*          JT has the following values and meanings: */
	/*           1 means a user-supplied full (NEQ by NEQ) Jacobian. */
	/*           2 means an internally generated (difference quotient) full */
	/*             Jacobian (using NEQ extra calls to F per df/dy value). */
	/*           4 means a user-supplied banded Jacobian. */
	/*           5 means an internally generated banded Jacobian (using */
	/*             ML+MU+1 extra calls to F per df/dy evaluation). */
	/*          If JT = 1 or 4, the user must supply a Subroutine JAC */
	/*          (the name is arbitrary) as described above under JAC. */
	/*          If JT = 2 or 5, a dummy argument can be used. */
	/* ----------------------------------------------------------------------- */
	/* Optional Inputs. */
	
	/* The following is a list of the optional inputs provided for in the */
	/* call sequence.  (See also Part 2.)  For each such input variable, */
	/* this table lists its name as used in this documentation, its */
	/* location in the call sequence, its meaning, and the default value. */
	/* The use of any of these inputs requires IOPT = 1, and in that */
	/* case all of these inputs are examined.  A value of zero for any */
	/* of these optional inputs will cause the default value to be used. */
	/* Thus to use a subset of the optional inputs, simply preload */
	/* locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and */
	/* then set those of interest to nonzero values. */
	
	/* Name    Location      Meaning and Default Value */
	
	/* H0      RWORK(5)  the step size to be attempted on the first step. */
	/*                   The default value is determined by the solver. */
	
	/* HMAX    RWORK(6)  the maximum absolute step size allowed. */
	/*                   The default value is infinite. */
	
	/* HMIN    RWORK(7)  the minimum absolute step size allowed. */
	/*                   The default value is 0.  (This lower bound is not */
	/*                   enforced on the final step before reaching TCRIT */
	/*                   when ITASK = 4 or 5.) */
	
	/* IXPR    IWORK(5)  flag to generate extra printing at method switches. */
	/*                   IXPR = 0 means no extra printing (the default). */
	/*                   IXPR = 1 means print data on each switch. */
	/*                   T, H, and NST will be printed on the same logical */
	/*                   unit as used for error messages. */
	
	/* MXSTEP  IWORK(6)  maximum number of (internally defined) steps */
	/*                   allowed during one call to the solver. */
	/*                   The default value is 500. */
	
	/* MXHNIL  IWORK(7)  maximum number of messages printed (per problem) */
	/*                   warning that T + H = T on a step (H = step size). */
	/*                   This must be positive to result in a non-default */
	/*                   value.  The default value is 10. */
	
	/* MXORDN  IWORK(8)  the maximum order to be allowed for the nonstiff */
	/*                   (Adams) method.  the default value is 12. */
	/*                   if MXORDN exceeds the default value, it will */
	/*                   be reduced to the default value. */
	/*                   MXORDN is held constant during the problem. */
	
	/* MXORDS  IWORK(9)  the maximum order to be allowed for the stiff */
	/*                   (BDF) method.  The default value is 5. */
	/*                   If MXORDS exceeds the default value, it will */
	/*                   be reduced to the default value. */
	/*                   MXORDS is held constant during the problem. */
	/* ----------------------------------------------------------------------- */
	/* Optional Outputs. */
	
	/* As optional additional output from DLSODA, the variables listed */
	/* below are quantities related to the performance of DLSODA */
	/* which are available to the user.  These are communicated by way of */
	/* the work arrays, but also have internal mnemonic names as shown. */
	/* except where stated otherwise, all of these outputs are defined */
	/* on any successful return from DLSODA, and on any return with */
	/* ISTATE = -1, -2, -4, -5, or -6.  On an illegal input return */
	/* (ISTATE = -3), they will be unchanged from their existing values */
	/* (if any), except possibly for TOLSF, LENRW, and LENIW. */
	/* On any error return, outputs relevant to the error will be defined, */
	/* as noted below. */
	
	/* Name    Location      Meaning */
	
	/* HU      RWORK(11) the step size in t last used (successfully). */
	
	/* HCUR    RWORK(12) the step size to be attempted on the next step. */
	
	/* TCUR    RWORK(13) the current value of the independent variable */
	/*                   which the solver has actually reached, i.e. the */
	/*                   current internal mesh point in t.  On output, TCUR */
	/*                   will always be at least as far as the argument */
	/*                   T, but may be farther (if interpolation was done). */
	
	/* TOLSF   RWORK(14) a tolerance scale factor, greater than 1.0, */
	/*                   computed when a request for too much accuracy was */
	/*                   detected (ISTATE = -3 if detected at the start of */
	/*                   the problem, ISTATE = -2 otherwise).  If ITOL is */
	/*                   left unaltered but RTOL and ATOL are uniformly */
	/*                   scaled up by a factor of TOLSF for the next call, */
	/*                   then the solver is deemed likely to succeed. */
	/*                   (The user may also ignore TOLSF and alter the */
	/*                   tolerance parameters in any other way appropriate.) */
	
	/* TSW     RWORK(15) the value of t at the time of the last method */
	/*                   switch, if any. */
	
	/* NST     IWORK(11) the number of steps taken for the problem so far. */
	
	/* NFE     IWORK(12) the number of f evaluations for the problem so far. */
	
	/* NJE     IWORK(13) the number of Jacobian evaluations (and of matrix */
	/*                   LU decompositions) for the problem so far. */
	
	/* NQU     IWORK(14) the method order last used (successfully). */
	
	/* NQCUR   IWORK(15) the order to be attempted on the next step. */
	
	/* IMXER   IWORK(16) the index of the component of largest magnitude in */
	/*                   the weighted local error vector ( E(i)/EWT(i) ), */
	/*                   on an error return with ISTATE = -4 or -5. */
	
	/* LENRW   IWORK(17) the length of RWORK actually required, assuming */
	/*                   that the length of RWORK is to be fixed for the */
	/*                   rest of the problem, and that switching may occur. */
	/*                   This is defined on normal returns and on an illegal */
	/*                   input return for insufficient storage. */
	
	/* LENIW   IWORK(18) the length of IWORK actually required, assuming */
	/*                   that the length of IWORK is to be fixed for the */
	/*                   rest of the problem, and that switching may occur. */
	/*                   This is defined on normal returns and on an illegal */
	/*                   input return for insufficient storage. */
	
	/* MUSED   IWORK(19) the method indicator for the last successful step: */
	/*                   1 means Adams (nonstiff), 2 means BDF (stiff). */
	
	/* MCUR    IWORK(20) the current method indicator: */
	/*                   1 means Adams (nonstiff), 2 means BDF (stiff). */
	/*                   This is the method to be attempted */
	/*                   on the next step.  Thus it differs from MUSED */
	/*                   only if a method switch has just been made. */
	
	/* The following two arrays are segments of the RWORK array which */
	/* may also be of interest to the user as optional outputs. */
	/* For each array, the table below gives its internal name, */
	/* its base address in RWORK, and its description. */
	
	/* Name    Base Address      Description */
	
	/* YH      21             the Nordsieck history array, of size NYH by */
	/*                        (NQCUR + 1), where NYH is the initial value */
	/*                        of NEQ.  For j = 0,1,...,NQCUR, column j+1 */
	/*                        of YH contains HCUR**j/factorial(j) times */
	/*                        the j-th derivative of the interpolating */
	/*                        polynomial currently representing the solution, */
	/*                        evaluated at T = TCUR. */
	
	/* ACOR     LACOR         array of size NEQ used for the accumulated */
	/*         (from Common   corrections on each step, scaled on output */
	/*           as noted)    to represent the estimated local error in y */
	/*                        on the last step.  This is the vector E in */
	/*                        the description of the error control.  It is */
	/*                        defined only on a successful return from */
	/*                        DLSODA.  The base address LACOR is obtained by */
	/*                        including in the user's program the */
	/*                        following 2 lines: */
	/*                           COMMON /DLS001/ RLS(218), ILS(37) */
	/*                           LACOR = ILS(22) */
	
	/* ----------------------------------------------------------------------- */
	/* Part 2.  Other Routines Callable. */
	
	/* The following are optional calls which the user may make to */
	/* gain additional capabilities in conjunction with DLSODA. */
	/* (The routines XSETUN and XSETF are designed to conform to the */
	/* SLATEC error handling package.) */
	
	/*     Form of Call                  Function */
	/*   CALL XSETUN(LUN)          set the logical unit number, LUN, for */
	/*                             output of messages from DLSODA, if */
	/*                             the default is not desired. */
	/*                             The default value of LUN is 6. */
	
	/*   CALL XSETF(MFLAG)         set a flag to control the printing of */
	/*                             messages by DLSODA. */
	/*                             MFLAG = 0 means do not print. (Danger: */
	/*                             This risks losing valuable information.) */
	/*                             MFLAG = 1 means print (the default). */
	
	/*                             Either of the above calls may be made at */
	/*                             any time and will take effect immediately. */
	
	/*   CALL DSRCMA(RSAV,ISAV,JOB) saves and restores the contents of */
	/*                             the internal Common blocks used by */
	/*                             DLSODA (see Part 3 below). */
	/*                             RSAV must be a real array of length 240 */
	/*                             or more, and ISAV must be an int */
	/*                             array of length 46 or more. */
	/*                             JOB=1 means save Common into RSAV/ISAV. */
	/*                             JOB=2 means restore Common from RSAV/ISAV. */
	/*                                DSRCMA is useful if one is */
	/*                             interrupting a run and restarting */
	/*                             later, or alternating between two or */
	/*                             more problems solved with DLSODA. */
	
	/*   CALL DINTDY(,,,,,)        provide derivatives of y, of various */
	/*        (see below)          orders, at a specified point t, if */
	/*                             desired.  It may be called only after */
	/*                             a successful return from DLSODA. */
	
	/* The detailed instructions for using DINTDY are as follows. */
	/* The form of the call is: */
	
	/*   CALL DINTDY (T, K, RWORK(21), NYH, DKY, IFLAG) */
	
	/* The input parameters are: */
	
	/* T         = value of independent variable where answers are desired */
	/*             (normally the same as the T last returned by DLSODA). */
	/*             For valid results, T must lie between TCUR - HU and TCUR. */
	/*             (See optional outputs for TCUR and HU.) */
	/* K         = int order of the derivative desired.  K must satisfy */
	/*             0 .le. K .le. NQCUR, where NQCUR is the current order */
	/*             (see optional outputs).  The capability corresponding */
	/*             to K = 0, i.e. computing y(T), is already provided */
	/*             by DLSODA directly.  Since NQCUR .ge. 1, the first */
	/*             derivative dy/dt is always available with DINTDY. */
	/* RWORK(21) = the base address of the history array YH. */
	/* NYH       = column length of YH, equal to the initial value of NEQ. */
	
	/* The output parameters are: */
	
	/* DKY       = a real array of length NEQ containing the computed value */
	/*             of the K-th derivative of y(t). */
	/* IFLAG     = int flag, returned as 0 if K and T were legal, */
	/*             -1 if K was illegal, and -2 if T was illegal. */
	/*             On an error return, a message is also written. */
	/* ----------------------------------------------------------------------- */
	/* Part 3.  Common Blocks. */
	
	/* If DLSODA is to be used in an overlay situation, the user */
	/* must declare, in the primary overlay, the variables in: */
	/*   (1) the call sequence to DLSODA, and */
	/*   (2) the two internal Common blocks */
	/*         /DLS001/  of length  255  (218 double  precision words */
	/*                      followed by 37 int words), */
	/*         /DLSA01/  of length  31    (22 double  precision words */
	/*                      followed by  9 int words). */
	
	/* If DLSODA is used on a system in which the contents of internal */
	/* Common blocks are not preserved between calls, the user should */
	/* declare the above Common blocks in the calling program to insure */
	/* that their contents are preserved. */
	
	/* If the solution of a given problem by DLSODA is to be interrupted */
	/* and then later continued, such as when restarting an interrupted run */
	/* or alternating between two or more problems, the user should save, */
	/* following the return from the last DLSODA call prior to the */
	/* interruption, the contents of the call sequence variables and the */
	/* internal Common blocks, and later restore these values before the */
	/* next DLSODA call for that problem.  To save and restore the Common */
	/* blocks, use Subroutine DSRCMA (see Part 2 above). */
	
	/* ----------------------------------------------------------------------- */
	/* Part 4.  Optionally Replaceable Solver Routines. */
	
	/* Below is a description of a routine in the DLSODA package which */
	/* relates to the measurement of errors, and can be */
	/* replaced by a user-supplied version, if desired.  However, since such */
	/* a replacement may have a major impact on performance, it should be */
	/* done only when absolutely necessary, and only with great caution. */
	/* (Note: The means by which the package version of a routine is */
	/* superseded by the user's version may be system-dependent.) */
	
	/* (a) DEWSET. */
	/* The following subroutine is called just before each internal */
	/* integration step, and sets the array of error weights, EWT, as */
	/* described under ITOL/RTOL/ATOL above: */
	/*     Subroutine DEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT) */
	/* where NEQ, ITOL, RTOL, and ATOL are as in the DLSODA call sequence, */
	/* YCUR contains the current dependent variable vector, and */
	/* EWT is the array of weights set by DEWSET. */
	
	/* If the user supplies this subroutine, it must return in EWT(i) */
	/* (i = 1,...,NEQ) a positive quantity suitable for comparing errors */
	/* in y(i) to.  The EWT array returned by DEWSET is passed to the */
	/* DMNORM routine, and also used by DLSODA in the computation */
	/* of the optional output IMXER, and the increments for difference */
	/* quotient Jacobians. */
	
	/* In the user-supplied version of DEWSET, it may be desirable to use */
	/* the current values of derivatives of y.  Derivatives up to order NQ */
	/* are available from the history array YH, described above under */
	/* optional outputs.  In DEWSET, YH is identical to the YCUR array, */
	/* extended to NQ + 1 columns with a column length of NYH and scale */
	/* factors of H**j/factorial(j).  On the first call for the problem, */
	/* given by NST = 0, NQ is 1 and H is temporarily set to 1.0. */
	/* NYH is the initial value of NEQ.  The quantities NQ, H, and NST */
	/* can be obtained by including in DEWSET the statements: */
	/*     DOUBLE PRECISION RLS */
	/*     COMMON /DLS001/ RLS(218),ILS(37) */
	/*     NQ = ILS(33) */
	/*     NST = ILS(34) */
	/*     H = RLS(212) */
	/* Thus, for example, the current value of dy/dt can be obtained as */
	/* YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is */
	/* unnecessary when NST = 0). */
	/* ----------------------------------------------------------------------- */
	
	/* ***REVISION HISTORY  (YYYYMMDD) */
	/* 19811102  DATE WRITTEN */
	/* 19820126  Fixed bug in tests of work space lengths; */
	/*           minor corrections in main prologue and comments. */
	/* 19870330  Major update: corrected comments throughout; */
	/*           removed TRET from Common; rewrote EWSET with 4 loops; */
	/*           fixed t test in INTDY; added Cray directives in STODA; */
	/*           in STODA, fixed DELP init. and logic around dprja_ call; */
	/*           combined routines to save/restore Common; */
	/*           passed LEVEL = 0 in error message calls (except run abort). */
	/* 19970225  Fixed lines setting JSTART = -2 in Subroutine LSODA. */
	/* 20010425  Major update: convert source lines to upper case; */
	/*           added *DECK lines; changed from 1 to * in dummy dimensions; */
	/*           changed names R1MACH/D1MACH to RUMACH/DUMACH; */
	/*           renamed routines for uniqueness across single/double  prec.; */
	/*           converted intrinsic names to generic form; */
	/*           removed ILLIN and NTREP (data loaded) from Common; */
	/*           removed all 'own' variables from Common; */
	/*           changed error messages to quoted strings; */
	/*           replaced XERRWV/XERRWD with 1993 revised version; */
	/*           converted prologues, comments, error messages to mixed case; */
	/*           numerous corrections to prologues and internal comments. */
	/* 20010507  Converted single precision source to double  precision. */
	/* 20010613  Revised excess accuracy test (to match rest of ODEPACK). */
	/* 20010808  Fixed bug in DPRJA (matrix in DBNORM call). */
	/* 20020502  Corrected declarations in descriptions of user routines. */
	/* 20031105  Restored 'own' variables to Common blocks, to enable */
	/*           interrupt/restart feature. */
	/* 20031112  Added SAVE statements for data-loaded constants. */
	
	/* ----------------------------------------------------------------------- */
	/* Other routines in the DLSODA package. */
	
	/* In addition to Subroutine DLSODA, the DLSODA package includes the */
	/* following subroutines and function routines: */
	/*  DINTDY   computes an interpolated value of the y vector at t = TOUT. */
	/*  DSTODA   is the core integrator, which does one step of the */
	/*           integration and the associated error control. */
	/*  DCFODE   sets all method coefficients and test constants. */
	/*  DPRJA    computes and preprocesses the Jacobian matrix J = df/dy */
	/*           and the Newton iteration matrix P = I - h*l0*J. */
	/*  DSOLSY   manages solution of linear system in chord iteration. */
	/*  DEWSET   sets the error weight vector EWT before each step. */
	/*  DMNORM   computes the weighted max-norm of a vector. */
	/*  DFNORM   computes the norm of a full matrix consistent with the */
	/*           weighted max-norm on vectors. */
	/*  DBNORM   computes the norm of a band matrix consistent with the */
	/*           weighted max-norm on vectors. */
	/*  DSRCMA   is a user-callable routine to save and restore */
	/*           the contents of the internal Common blocks. */
	/*  DGEFA and DGESL   are routines from LINPACK for solving full */
	/*           systems of linear algebraic equations. */
	/*  DGBFA and DGBSL   are routines from LINPACK for solving banded */
	/*           linear systems. */
	/*  DUMACH   computes the unit roundoff in a machine-independent manner. */
	/*  XERRWD, XSETUN, XSETF, IXSAV, and IUMACH  handle the printing of all */
	/*           error messages and warnings.  XERRWD is machine-dependent. */
	/* Note:  DMNORM, DFNORM, DBNORM, DUMACH, IXSAV, and IUMACH are */
	/* function routines.  All the others are subroutines. */
	
	/* ----------------------------------------------------------------------- */
	/* ----------------------------------------------------------------------- */
	/* The following two internal Common blocks contain */
	/* (a) variables which are local to any subroutine but whose values must */
	/*     be preserved between calls to the routine ("own" variables), and */
	/* (b) variables which are communicated between subroutines. */
	/* The block DLS001 is declared in subroutines DLSODA, DINTDY, DSTODA, */
	/* DPRJA, and DSOLSY. */
	/* The block DLSA01 is declared in subroutines DLSODA, DSTODA, and DPRJA. */
	/* Groups of variables are replaced by dummy arrays in the Common */
	/* declarations in routines where those variables are not used. */
	/* ----------------------------------------------------------------------- */
	
	
    /* Parameter adjustments */
    //--neq;
    //--y;
    //--rtol;
    //--atol;
    //--rwork;
   // --iwork;
	
    /* Function Body */
	/* ----------------------------------------------------------------------- */
	/* Block A. */
	/* This code block is executed on every call. */
	/* It tests ISTATE and ITASK for legality and branches appropriately. */
	/* If ISTATE .gt. 1 but the flag INIT shows that initialization has */
	/* not yet been done, an error return occurs. */
	/* If ISTATE = 1 and TOUT = T, return immediately. */
	/* ----------------------------------------------------------------------- */
    if (*istate < 1 || *istate > 3) {
		goto L601;
    }
    if (*itask < 1 || *itask > 5) {
		goto L602;
    }
    if (*istate == 1) {
		goto L10;
    }
    if (init == 0) {
		goto L603;
    }
    if (*istate == 2) {
		goto L200;
    }
    goto L20;
L10:
    init = 0;
    if (*tout == *t) {
		return errorCode;
    }
	/* ----------------------------------------------------------------------- */
	/* Block B. */
	/* The next code block is executed for the initial call (ISTATE = 1), */
	/* or for a continuation call with parameter changes (ISTATE = 3). */
	/* It contains checking of all inputs and various initializations. */
	
	/* First check legality of the non-optional inputs NEQ, ITOL, IOPT, */
	/* JT, ML, and MU. */
	/* ----------------------------------------------------------------------- */
L20:
    if (neq[0] <= 0) { //fixed
		goto L604;
    }
    if (*istate == 1) {
		goto L25;
    }
    if (neq[0] > n) { //fixed
		goto L605;
    }
L25:
    n = neq[0];  //fixed
    if (*itol < 1 || *itol > 4) {
		goto L606;
    }
    if (*iopt < 0 || *iopt > 1) {
		goto L607;
    }
    if (*jt == 3 || *jt < 1 || *jt > 5) {
		goto L608;
    }
    jtyp = *jt;
    if (*jt <= 2) {
		goto L30;
    }
    ml = iwork[0];
    mu = iwork[1];
    if (ml < 0 || ml >= n) {
		goto L609;
    }
    if (mu < 0 || mu >= n) {
		goto L610;
    }
L30:
	/* Next process and check the optional inputs. -------------------------- */
    if (*iopt == 1) {
		goto L40;
    }
    ixpr = 0;
    mxstep = mxstp0;
    mxhnil = mxhnl0;
    hmxi = 0.;
    hmin = 0.;
    if (*istate != 1) {
		goto L60;
    }
    h0 = 0.;
    mxordn = mord[0];
    mxords = mord[1];
    goto L60;
L40:
    ixpr = iwork[4];
    if (ixpr < 0 || ixpr > 1) {
		goto L611;
    }
    mxstep = iwork[5];
    if (mxstep < 0) {
		goto L612;
    }
    if (mxstep == 0) {
		mxstep = mxstp0;
    }
    mxhnil = iwork[6];
    if (mxhnil < 0) {
		goto L613;
    }
    if (mxhnil == 0) {
		mxhnil = mxhnl0;
    }
    if (*istate != 1) {
		goto L50;
    }
    h0 = rwork[4];
    mxordn = iwork[7];
    if (mxordn < 0) {
		goto L628;
    }
    if (mxordn == 0) {
		mxordn = 100;
    }
    mxordn = min(mxordn,mord[0]);
    mxords = iwork[8];
    if (mxords < 0) {
		goto L629;
    }
    if (mxords == 0) {
		mxords = 100;
    }
    mxords = min(mxords,mord[1]);
    if ((*tout - *t) * h0 < 0.) {
		goto L614;
    }
L50:
    hmax = rwork[5];
    if (hmax < 0.) {
		goto L615;
    }
    hmxi = 0.;
    if (hmax > 0.) {
		hmxi = 1. / hmax;
    }
    hmin = rwork[6];
    if (hmin < 0.) {
		goto L616;
    }
	/* ----------------------------------------------------------------------- */
	/* Set work array pointers and check lengths LRW and LIW. */
	/* If ISTATE = 1, METH is initialized to 1 here to facilitate the */
	/* checking of work space lengths. */
	/* Pointers to segments of RWORK and IWORK are named by prefixing L to */
	/* the name of the segment.  E.g., the segment YH starts at RWORK(LYH). */
	/* Segments of RWORK (in order) are denoted  YH, WM, EWT, SAVF, ACOR. */
	/* If the lengths provided are insufficient for the current method, */
	/* an error return occurs.  This is treated as illegal input on the */
	/* first call, but as a problem interruption with ISTATE = -7 on a */
	/* continuation call.  If the lengths are sufficient for the current */
	/* method but not for both methods, a warning message is sent. */
	/* ----------------------------------------------------------------------- */
L60:
    if (*istate == 1) {
		meth = 1;
    }
    if (*istate == 1) {
		nyh = n;
    }
    lyh = 21;
    len1n = (mxordn + 1) * nyh + 20;
    len1s = (mxords + 1) * nyh + 20;
    lwm = len1s + 1;
    if (*jt <= 2) {
		lenwm = n * n + 2;
    }
    if (*jt >= 4) {
		lenwm = ((ml << 1) + mu + 1) * n + 2;
    }
    len1s += lenwm;
    len1c = len1n;
    if (meth == 2) {
		len1c = len1s;
    }
    len1 = max(len1n,len1s);
    len2 = n * 3;
    lenrw = len1 + len2;
    lenrwc = len1c + len2;
    iwork[16] = lenrw;
    liwm = 1;
    leniw = n + 20;
    leniwc = 20;
    if (meth == 2) {
		leniwc = leniw;
    }
    iwork[17] = leniw;
    if (*istate == 1 && *lrw < lenrwc) {
		goto L617;
    }
    if (*istate == 1 && *liw < leniwc) {
		goto L618;
    }
    if (*istate == 3 && *lrw < lenrwc) {
		goto L550;
    }
    if (*istate == 3 && *liw < leniwc) {
		goto L555;
    }
    lewt = len1 + 1;
    insufr = 0;
    if (*lrw >= lenrw) {
		goto L65;
    }
    insufr = 2;
    lewt = len1c + 1;

#if (GMS_CUDA_DEBUG_ON) == 1
	 fprintf(stderr, "DLSODA-  Warning.. RWORK length is sufficient for now, but  \n");

	 fprintf(stderr, "      may not be later.  Integration will proceed anyway.   \n");

	 fprintf(stderr, "      Length needed is LENRW = I1, while LRW = I2.\n");
#endif
	 errorCode = 60;
	 
L65:
    lsavf = lewt + n;
    lacor = lsavf + n;
    insufi = 0;
    if (*liw >= leniw) {
		goto L70;
    }
    insufi = 2;
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  Warning.. IWORK length is sufficient for now, but  \n");

    fprintf(stderr, "      may not be later.  Integration will proceed anyway.   \n");

    fprintf(stderr, "      Length needed is LENIW = I1, while LIW = I2.\n");
#endif
	errorCode = 65;
	
L70:
	/* Check RTOL and ATOL for legality. ------------------------------------ */
    rtoli = rtol[0];
    atoli = atol[0];
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		if (*itol >= 3) {
			rtoli = rtol[i__ -1];
		}
		if (*itol == 2 || *itol == 4) {
			atoli = atol[i__ -1];
		}
		if (rtoli < 0.) {
			goto L619;
		}
		if (atoli < 0.) {
			goto L620;
		}
		/* L75: */
    }
    if (*istate == 1) {
		goto L100;
    }
	/* If ISTATE = 3, set flag to signal parameter changes to DSTODA. ------- */
    jstart = -1;
    if (n == nyh) {
		goto L200;
    }
	/* NEQ was reduced.  Zero part of YH to avoid undefined references. ----- */
    i1 = lyh + l * nyh;
    i2 = lyh + (maxord + 1) * nyh - 1;
    if (i1 > i2) {
		goto L200;
    }
    i__1 = i2;
    for (i__ = i1; i__ <= i__1; ++i__) {
		/* L95: */
		rwork[i__ -1] = 0.;
    }
    goto L200;
	/* ----------------------------------------------------------------------- */
	/* Block C. */
	/* The next block is for the initial call only (ISTATE = 1). */
	/* It contains all remaining initializations, the initial call to F, */
	/* and the calculation of the initial step size. */
	/* The error weights in EWT are inverted after being loaded. */
	/* ----------------------------------------------------------------------- */
L100:
    uround = dumach_(common);
    tn = *t;
    tsw = *t;
    maxord = mxordn;
    if (*itask != 4 && *itask != 5) {
		goto L110;
    }
    tcrit = rwork[0];
    if ((tcrit - *tout) * (*tout - *t) < 0.) {
		goto L625;
    }
    if (h0 != 0. && (*t + h0 - tcrit) * h0 > 0.) {
		h0 = tcrit - *t;
    }
L110:
    jstart = 0;
    nhnil = 0;
    nst = 0;
    nje = 0;
    nslast = 0;
    hu = 0.;
    nqu = 0;
    mused = 0;
    miter = 0;
    ccmax = .3;
    maxcor = 3;
    msbp = 20;
    mxncf = 10;
	/* Initial call to F.  (LF0 points to YH(*,2).) ------------------------- */
    lf0 = lyh + nyh;
    f(neq, t, y, &rwork[lf0 -1]); //fixed neq y
    nfe = 1;
	/* Load the initial value vector in YH. --------------------------------- */
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L115: */
		rwork[i__ + lyh - 1 -1] = y[i__ - 1];
    }
	/* Load and invert the EWT array.  (H is temporarily set to 1.0.) ------- */
    nq = 1;
    h__ = 1.;
    dewset_(&n, itol, rtol, atol, &rwork[lyh -1], &rwork[lewt -1], common);
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		if (rwork[i__ + lewt - 1 -1] <= 0.) {
			goto L621;
		}
		/* L120: */
		rwork[i__ + lewt - 1 -1] = 1. / rwork[i__ + lewt - 1 -1];
    }
	/* ----------------------------------------------------------------------- */
	/* The coding below computes the step size, H0, to be attempted on the */
	/* first step, unless the user has supplied a value for this. */
	/* First check that TOUT - T differs significantly from zero. */
	/* A scalar tolerance quantity TOL is computed, as MAX(RTOL(i)) */
	/* if this is positive, or MAX(ATOL(i)/ABS(Y(i))) otherwise, adjusted */
	/* so as to be between 100*UROUND and 1.0E-3. */
	/* Then the computed value H0 is given by: */
	
	/*   H0**(-2)  =  1./(TOL * w0**2)  +  TOL * (norm(F))**2 */
	
	/* where   w0     = MAX ( ABS(T), ABS(TOUT) ), */
	/*         F      = the initial value of the vector f(t,y), and */
	/*         norm() = the weighted vector norm used throughout, given by */
	/*                  the DMNORM function routine, and weighted by the */
	/*                  tolerances initially loaded into the EWT array. */
	/* The sign of H0 is inferred from the initial values of TOUT and T. */
	/* ABS(H0) is made .le. ABS(TOUT-T) in any case. */
	/* ----------------------------------------------------------------------- */
    if (h0 != 0.) {
		goto L180;
    }
    tdist = (d__1 = *tout - *t, fabs(d__1));
	/* Computing MAX */
    d__1 = fabs(*t), d__2 = fabs(*tout);
    w0 = max(d__1,d__2);
    if (tdist < uround * 2. * w0) {
		goto L622;
    }
    tol = rtol[0];
    if (*itol <= 2) {
		goto L140;
    }
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L130: */
		/* Computing MAX */
		d__1 = tol, d__2 = rtol[i__ -1];
		tol = max(d__1,d__2);
    }
L140:
    if (tol > 0.) {
		goto L160;
    }
    atoli = atol[0];
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		if (*itol == 2 || *itol == 4) {
			atoli = atol[i__ -1];
		}
		ayi = (d__1 = y[i__ - 1], fabs(d__1));
		if (ayi != 0.) {
			/* Computing MAX */
			d__1 = tol, d__2 = atoli / ayi;
			tol = max(d__1,d__2);
		}
		/* L150: */
    }
L160:
	/* Computing MAX */
    d__1 = tol, d__2 = uround * 100.;
    tol = max(d__1,d__2);
    tol = min(tol,.001);
    sum = dmnorm_(&n, &rwork[lf0 -1], &rwork[lewt -1], common);
	/* Computing 2nd power */
    d__1 = sum;
    sum = 1. / (tol * w0 * w0) + tol * (d__1 * d__1);
    h0 = 1. / sqrt(sum);
    h0 = min(h0,tdist);
    d__1 = *tout - *t;
    h0 = d_sign(&h0, &d__1);
	/* Adjust H0 if necessary to meet HMAX bound. --------------------------- */
L180:
    rh = fabs(h0) * hmxi;
    if (rh > 1.) {
		h0 /= rh;
    }
	/* Load H with H0 and scale YH(*,2) by H0. ------------------------------ */
    h__ = h0;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L190: */
		rwork[i__ + lf0 - 1 -1] = h0 * rwork[i__ + lf0 - 1 -1];
    }
    goto L270;
	/* ----------------------------------------------------------------------- */
	/* Block D. */
	/* The next code block is for continuation calls only (ISTATE = 2 or 3) */
	/* and is to check stop conditions before taking a step. */
	/* ----------------------------------------------------------------------- */
L200:
    nslast = nst;
    switch (*itask) {
		case 1:  goto L210;
		case 2:  goto L250;
		case 3:  goto L220;
		case 4:  goto L230;
		case 5:  goto L240;
    }
L210:
    if ((tn - *tout) * h__ < 0.) {
		goto L250;
    }
    dintdy_(tout, 0, &rwork[lyh -1], &nyh, y, &iflag, common); //fixed y
    if (iflag != 0) {
		goto L627;
    }
    *t = *tout;
    goto L420;
L220:
    tp = tn - hu * (uround * 100. + 1.);
    if ((tp - *tout) * h__ > 0.) {
		goto L623;
    }
    if ((tn - *tout) * h__ < 0.) {
		goto L250;
    }
    *t = tn;
    goto L400;
L230:
    tcrit = rwork[0];
    if ((tn - tcrit) * h__ > 0.) {
		goto L624;
    }
    if ((tcrit - *tout) * h__ < 0.) {
		goto L625;
    }
    if ((tn - *tout) * h__ < 0.) {
		goto L245;
    }
    dintdy_(tout, 0, &rwork[lyh -1], &nyh, y, &iflag, common); //fixed y
    if (iflag != 0) {
		goto L627;
    }
    *t = *tout;
    goto L420;
L240:
    tcrit = rwork[0];
    if ((tn - tcrit) * h__ > 0.) {
		goto L624;
    }
L245:
    hmx = fabs(tn) + fabs(h__);
    ihit = (d__1 = tn - tcrit, fabs(d__1)) <= uround * 100. *
	hmx;
    if (ihit) {
		*t = tcrit;
    }
    if (ihit) {
		goto L400;
    }
    tnext = tn + h__ * (uround * 4. + 1.);
    if ((tnext - tcrit) * h__ <= 0.) {
		goto L250;
    }
    h__ = (tcrit - tn) * (1. - uround * 4.);
    if (*istate == 2 && jstart >= 0) {
		jstart = -2;
    }
	/* ----------------------------------------------------------------------- */
	/* Block E. */
	/* The next block is normally executed for all calls and contains */
	/* the call to the one-step core integrator DSTODA. */
	
	/* This is a looping point for the integration steps. */
	
	/* First check for too many steps being taken, update EWT (if not at */
	/* start of problem), check for too much accuracy being requested, and */
	/* check for H below the roundoff level in T. */
	/* ----------------------------------------------------------------------- */
L250:
    if (meth == mused) {
		goto L255;
    }
    if (insufr == 1) {
		goto L550;
    }
    if (insufi == 1) {
		goto L555;
    }
L255:
    if (nst - nslast >= mxstep) {
		goto L500;
    }
    dewset_(&n, itol, rtol, atol, &rwork[lyh -1], &rwork[lewt -1], common);
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		if (rwork[i__ + lewt - 1 -1] <= 0.) {
			goto L510;
		}
		/* L260: */
		rwork[i__ + lewt - 1 -1] = 1. / rwork[i__ + lewt - 1 -1];
    }
L270:
    tolsf = uround * dmnorm_(&n, &rwork[lyh -1], &rwork[lewt -1], common);
    if (tolsf <= 1.) {
		goto L280;
    }
    tolsf *= 2.;
    if (nst == 0) {
		goto L626;
    }
    goto L520;
L280:
    if (tn + h__ != tn) {
		goto L290;
    }
    ++nhnil;
    if (nhnil > mxhnil) {
		goto L290;
    }
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  Warning..Internal T (=R1) and H (=R2) are\n");

    fprintf(stderr, "      such that in the machine, T + H = T on the next step  \n");

    fprintf(stderr, "     (H = step size). Solver will continue anyway.\n");
#endif
	errorCode = 280;
    if (nhnil < mxhnil) {
		goto L290;
    }
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  Above warning has been issued I1 times.  \n");

    fprintf(stderr, "     It will not be issued again for this problem.\n");
#endif
errorCode = 285;
   
L290:
	/* ----------------------------------------------------------------------- */
	/*   CALL DSTODA(NEQ,Y,YH,NYH,YH,EWT,SAVF,ACOR,WM,IWM,F,JAC,DPRJA,DSOLSY) */
	/* ----------------------------------------------------------------------- */
    dstoda_(neq, y, &rwork[lyh -1], &nyh, &rwork[lyh -1], &rwork[lewt -1], &rwork[lsavf -1], &rwork[lacor -1], &rwork[lwm -1], &iwork[liwm -1], f, jac, common); //fixed neq y
    kgo = 1 - kflag;
    switch (kgo) {
		case 1:  goto L300;
		case 2:  goto L530;
		case 3:  goto L540;
    }
	/* ----------------------------------------------------------------------- */
	/* Block F. */
	/* The following block handles the case of a successful return from the */
	/* core integrator (KFLAG = 0). */
	/* If a method switch was just made, record TSW, reset MAXORD, */
	/* set JSTART to -1 to signal DSTODA to complete the switch, */
	/* and do extra printing of data if IXPR = 1. */
	/* Then, in any case, check for stop conditions. */
	/* ----------------------------------------------------------------------- */
L300:
    init = 1;
    if (meth == mused) {
		goto L310;
    }
    tsw = tn;
    maxord = mxordn;
    if (meth == 2) {
		maxord = mxords;
    }
    if (meth == 2) {
		rwork[lwm -1] = sqrt(uround);
    }
    insufr = min(insufr,1);
    insufi = min(insufi,1);
    jstart = -1;
    if (ixpr == 0) {
		goto L310;
    }
    if (meth == 2) {

		fprintf(stderr, "DLSODA- A switch to the BDF (stiff) method has occurred\n");

		//xerrwd_(msg, &c__60, 105, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b62, &c_b62, (ftnlen)60);
    }
    if (meth == 1) {
#if (GMS_CUDA_DEBUG_ON) == 1
		fprintf(stderr, "DLSODA- A switch to the Adams (nonstiff) method has occ\n");
#endif
		//xerrwd_(msg, &c__60, 106, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b62, &c_b62, (ftnlen)60);
    }
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "     at T = R1,  tentative step size H = R2,  step NST = I1 \n");
#endif
L310:
    switch (*itask) {
		case 1:  goto L320;
		case 2:  goto L400;
		case 3:  goto L330;
		case 4:  goto L340;
		case 5:  goto L350;
    }
	/* ITASK = 1.  If TOUT has been reached, interpolate. ------------------- */
L320:
    if ((tn - *tout) * h__ < 0.) {
		goto L250;
    }
    dintdy_(tout, 0, &rwork[lyh -1], &nyh, y, &iflag, common);
    *t = *tout;
    goto L420;
	/* ITASK = 3.  Jump to exit if TOUT was reached. ------------------------ */
L330:
    if ((tn - *tout) * h__ >= 0.) {
		goto L400;
    }
    goto L250;
	/* ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary. */
L340:
    if ((tn - *tout) * h__ < 0.) {
		goto L345;
    }
    dintdy_(tout, 0, &rwork[lyh -1], &nyh, y, &iflag, common);
    *t = *tout;
    goto L420;
L345:
    hmx = fabs(tn) + fabs(h__);
    ihit = (d__1 = tn - tcrit, fabs(d__1)) <= uround * 100. *
	hmx;
    if (ihit) {
		goto L400;
    }
    tnext = tn + h__ * (uround * 4. + 1.);
    if ((tnext - tcrit) * h__ <= 0.) {
		goto L250;
    }
    h__ = (tcrit - tn) * (1. - uround * 4.);
    if (jstart >= 0) {
		jstart = -2;
    }
    goto L250;
	/* ITASK = 5.  See if TCRIT was reached and jump to exit. --------------- */
L350:
    hmx = fabs(tn) + fabs(h__);
    ihit = (d__1 = tn - tcrit, fabs(d__1)) <= uround * 100. *
	hmx;
	/* ----------------------------------------------------------------------- */
	/* Block G. */
	/* The following block handles all successful returns from DLSODA. */
	/* If ITASK .ne. 1, Y is loaded from YH and T is set accordingly. */
	/* ISTATE is set to 2, and the optional outputs are loaded into the */
	/* work arrays before returning. */
	/* ----------------------------------------------------------------------- */
L400:
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L410: */
		y[i__ - 1] = rwork[i__ + lyh - 1 -1];  //fixed y
    }
    *t = tn;
    if (*itask != 4 && *itask != 5) {
		goto L420;
    }
    if (ihit) {
		*t = tcrit;
    }
L420:
    *istate = 2;
    rwork[10] = hu;
    rwork[11] = h__;
    rwork[12] = tn;
    rwork[14] = tsw;
    iwork[10] = nst;
    iwork[11] = nfe;
    iwork[12] = nje;
    iwork[13] = nqu;
    iwork[14] = nq;
    iwork[18] = mused;
    iwork[19] = meth;
	
	return errorCode;
	/* ----------------------------------------------------------------------- */
	/* Block H. */
	/* The following block handles all unsuccessful returns other than */
	/* those for illegal input.  First the error message routine is called. */
	/* If there was an error test or convergence test failure, IMXER is set. */
	/* Then Y is loaded from YH and T is set to TN. */
	/* The optional outputs are loaded into the work arrays before returning. */
	/* ----------------------------------------------------------------------- */
	/* The maximum number of steps was taken before reaching TOUT. ---------- */
L500:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  At current T (=R1), MXSTEP (=I1) steps   \n");

    fprintf(stderr, "      taken on this call before reaching TOUT     \n");
#endif
	errorCode = 500;
    *istate = -1;
    goto L580;
	/* EWT(i) .le. 0.0 for some i (not at start of problem). ---------------- */
L510:
    ewti = rwork[lewt + i__ - 1 -1];
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  At T (=R1), EWT(%d) has become R2 .le. 0.\n",ewti);
#endif
    errorCode = 510;
	*istate = -6;
    goto L580;
	/* Too much accuracy requested for machine precision. ------------------- */
L520:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  At T (=R1), too much accuracy requested  \n");

    fprintf(stderr, "      for precision of machine..  See TOLSF (=R2) \n");
#endif
	errorCode = 520;
    rwork[13] = tolsf;
    *istate = -2;
    goto L580;
	/* KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. ----- */
L530:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  At T(=R1) and step size H(=R2), the error\n");

    fprintf(stderr, "      test failed repeatedly or with ABS(H) = HMIN\n");
#endif
    
    errorCode = 530;
	*istate = -4;
    goto L560;
	/* KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ---- */
L540:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  At T (=R1) and step size H (=R2), the    \n");

    fprintf(stderr, "      corrector convergence failed repeatedly     \n");

    fprintf(stderr, "      or with ABS(H) = HMIN   \n");
#endif
	errorCode = 540;
    *istate = -5;
    goto L560;
	/* RWORK length too small to proceed. ----------------------------------- */
L550:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  At current T(=R1), RWORK length too small\n");

    fprintf(stderr, "      to proceed.  The integration was otherwise successful.\n");
#endif
    errorCode = 550;
	*istate = -7;
    goto L580;
	/* IWORK length too small to proceed. ----------------------------------- */
L555:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  At current T(=R1), IWORK length too small\n");

    fprintf(stderr, "      to proceed.  The integration was otherwise successful.\n");
#endif
   
    errorCode = 555;
	*istate = -7;
    goto L580;
	/* Compute IMXER if relevant. ------------------------------------------- */
L560:
    big = 0.;
    imxer = 1;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		size = (d__1 = rwork[i__ + lacor - 1 -1] * rwork[i__ + lewt - 1 -1], fabs(d__1));
		if (big >= size) {
			goto L570;
		}
		big = size;
		imxer = i__;
	L570:
		;
    }
    iwork[15] = imxer;
	/* Set Y vector, T, and optional outputs. ------------------------------- */
L580:
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L590: */
		y[i__ - 1] = rwork[i__ + lyh - 1 -1]; //fixed y
    }
    *t = tn;
    rwork[10] = hu;
    rwork[11] = h__;
    rwork[12] = tn;
    rwork[14] = tsw;
    iwork[10] = nst;
    iwork[11] = nfe;
    iwork[12] = nje;
    iwork[13] = nqu;
    iwork[14] = nq;
    iwork[18] = mused;
    iwork[19] = meth;
    
	return errorCode;
	/* ----------------------------------------------------------------------- */
	/* Block I. */
	/* The following block handles all error returns due to illegal input */
	/* (ISTATE = -3), as detected before calling the core integrator. */
	/* First the error message routine is called.  If the illegal input */
	/* is a negative ISTATE, the run is aborted (apparent infinite loop). */
	/* ----------------------------------------------------------------------- */
L601:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  ISTATE (=I1) illegal.\n");
#endif
    if (*istate < 0) {
		errorCode = 6011;
		goto L800;
    }
    errorCode = 601;
	goto L700;
L602:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  ITASK (=I1) illegal. \n");
#endif
    errorCode = 602;
	goto L700;
L603:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  ISTATE .gt. 1 but DLSODA not initialized.\n");
#endif
    errorCode = 603;
	goto L700;
L604:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  NEQ (=I1) .lt. 1     \n");
#endif
    errorCode = 604;
	goto L700;
L605:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  ISTATE = 3 and NEQ increased (I1 to I2). \n");
#endif
    errorCode = 605;
	goto L700;
L606:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  ITOL (=I1) illegal.  \n");
#endif
    errorCode = 606;
	goto L700;
L607:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  IOPT (=I1) illegal.  \n");
#endif
    errorCode = 607;
	goto L700;
L608:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  JT (=I1) illegal.    \n");
#endif
    errorCode = 608;
	goto L700;
L609:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  ML (=I1) illegal: .lt.0 or .ge.NEQ (=I2) \n");
#endif
    errorCode = 609;
	goto L700;
L610:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  MU (=I1) illegal: .lt.0 or .ge.NEQ (=I2) \n");
#endif
   errorCode = 610;
	goto L700;
L611:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  IXPR (=I1) illegal.  \n");
#endif
    errorCode = 611;
	goto L700;
L612:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  MXSTEP (=I1) .lt. 0  \n");
#endif
    errorCode = 612;
	goto L700;
L613:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  MXHNIL (=I1) .lt. 0  \n");
#endif
	errorCode = 613;
    goto L700;
L614:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  TOUT (=R1) behind T (=R2)      \n");

    fprintf(stderr, "      Integration direction is given by H0 (=R1)  \n");
#endif
    errorCode = 614;
	goto L700;
L615:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  HMAX (=R1) .lt. 0.0  \n");
#endif
    errorCode = 615;
	goto L700;
L616:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  HMIN (=R1) .lt. 0.0  \n");
#endif
	goto L700;
L617:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)\n");
#endif
    errorCode = 617;
	goto L700;
L618:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)\n");
#endif
    errorCode = 618;
	goto L700;
L619:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  RTOL(I1) is R1 .lt. 0.0        \n");
#endif
    errorCode = 619;
	goto L700;
L620:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  ATOL(I1) is R1 .lt. 0.0        \n");
#endif
    errorCode = 620;
	goto L700;
L621:
    ewti = rwork[lewt + i__ - 1 -1];
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  EWT(I1) is R1 .le. 0.0         \n");
#endif
    errorCode = 621;
	goto L700;
L622:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  TOUT(=R1) too close to T(=R2) to start integration.\n");
#endif
    errorCode = 622;
	goto L700;
L623:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  \n");
#endif
    errorCode = 623;
	goto L700;
L624:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   \n");
#endif
    errorCode = 624;
	goto L700;
L625:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   \n");
#endif
    errorCode = 625;
	goto L700;
L626:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  At start of problem, too much accuracy   \n");

    fprintf(stderr, "      requested for precision of machine..  See TOLSF (=R1) \n");
#endif
    errorCode = 626;
	rwork[13] = tolsf;
    goto L700;
L627:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  Trouble in DINTDY.  ITASK = I1, TOUT = R1\n");
#endif
    errorCode = 627;
	goto L700;
L628:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  MXORDN (=I1) .lt. 0  \n");
#endif
    errorCode = 628;
	goto L700;
L629:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  MXORDS (=I1) .lt. 0  \n");
#endif
L700:
    *istate = -3;
	return errorCode;
	
L800:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DLSODA-  Run aborted.. apparent infinite loop.    \n");
#endif
	return errorCode;
	/* ----------------------- End of Subroutine DLSODA ---------------------- */
} /* dlsoda_ */

/* DECK DINTDY */
/* Subroutine */ 
__device__ int dintdy_(double  * __restrict__ t, 
                       int k, 
                       double  * __restrict__ yh, 
                       int * __restrict__ NOT_nyh, 
                       double  * __restrict__ dky, 
                       int * __restrict__ iflag, 
                       struct cuLsodaCommonBlock *__restrict__ common)
{
    /* System generated locals */
    int yh_dim1 = 0;
	int yh_offset = 0;
	int i__1 = 0;
	int i__2 = 0;
    double  d__1 = 0.;
	
    /* Builtin functions */

	
	
	
	
    /* Local variables */
     double  c__ = 0.;
     int i__ = 0;
	 int j = 0;
     double  r__ = 0;
	 double  s = 0; 
     int ic = 0;
	 int jb = 0;
	 int jj = 0;
     double  tp = 0;
     int jb2 = 0;
	 int jj1 = 0;
	 int jp1 = 0;
	
	/* ***BEGIN PROLOGUE  DINTDY */
	/* ***SUBSIDIARY */
	/* ***PURPOSE  Interpolate solution derivatives. */
	/* ***TYPE      DOUBLE PRECISION (SINTDY-S, DINTDY-D) */
	/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
	/* ***DESCRIPTION */
	
	/*  DINTDY computes interpolated values of the K-th derivative of the */
	/*  dependent variable vector y, and stores it in DKY.  This routine */
	/*  is called within the package with K = 0 and T = TOUT, but may */
	/*  also be called by the user for any K up to the current order. */
	/*  (See detailed instructions in the usage documentation.) */
	
	/*  The computed values in DKY are gotten by interpolation using the */
	/*  Nordsieck history array YH.  This array corresponds uniquely to a */
	/*  vector-valued polynomial of degree NQCUR or less, and DKY is set */
	/*  to the K-th derivative of this polynomial at T. */
	/*  The formula for DKY is: */
	/*               q */
	/*   DKY(i)  =  sum  c(j,K) * (T - tn)**(j-K) * h**(-j) * YH(i,j+1) */
	/*              j=K */
	/*  where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, tn = TCUR, h = HCUR. */
	/*  The quantities  nq = NQCUR, l = nq+1, N = NEQ, tn, and h are */
	/*  communicated by COMMON.  The above sum is done in reverse order. */
	/*  IFLAG is returned negative if either K or T is out of bounds. */
	
	/* ***SEE ALSO  DLSODE */
	/* ***ROUTINES CALLED  XERRWD */
	/* ***COMMON BLOCKS    DLS001 */
	/* ***REVISION HISTORY  (YYMMDD) */
	/*   791129  DATE WRITTEN */
	/*   890501  Modified prologue to SLATEC/LDOC format.  (FNF) */
	/*   890503  Minor cosmetic changes.  (FNF) */
	/*   930809  Renamed to allow single/double  precision versions. (ACH) */
	/*   010418  Reduced size of Common block /DLS001/. (ACH) */
	/*   031105  Restored 'own' variables to Common block /DLS001/, to */
	/*           enable interrupt/restart feature. (ACH) */
	/*   050427  Corrected roundoff decrement in TP. (ACH) */
	/* ***END PROLOGUE  DINTDY */
	/* **End */
	
	/* ***FIRST EXECUTABLE STATEMENT  DINTDY */
    /* Parameter adjustments */
    yh_dim1 = *NOT_nyh;
    yh_offset = 1 + yh_dim1;
    //yh -= yh_offset;
    //--dky;
	
    /* Function Body */
    *iflag = 0;
    if (k < 0 || k > nq) {
		goto L80;
    }
    d__1 = fabs(tn) + fabs(hu);
    tp = tn - hu - uround * 100. * d_sign(&d__1, &hu);
    if ((*t - tp) * (*t - tn) > 0.) {
		goto L90;
    }
	
    s = (*t - tn) / h__;
    ic = 1;
    if (k == 0) {
		goto L15;
    }
    jj1 = l - k;
    i__1 = nq;
    for (jj = jj1; jj <= i__1; ++jj) {
		/* L10: */
		ic *= jj;
    }
L15:
    c__ = (double ) ic;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L20: */
		dky[i__ -1] = c__ * yh[i__ + l * yh_dim1 -yh_offset];
    }
    if (k == nq) {
		goto L55;
    }
    jb2 = nq - k;
    i__1 = jb2;
    for (jb = 1; jb <= i__1; ++jb) {
		j = nq - jb;
		jp1 = j + 1;
		ic = 1;
		if (k == 0) {
			goto L35;
		}
		jj1 = jp1 - k;
		i__2 = j;
		for (jj = jj1; jj <= i__2; ++jj) {
			/* L30: */
			ic *= jj;
		}
	L35:
		c__ = (double ) ic;
		i__2 = n;
		for (i__ = 1; i__ <= i__2; ++i__) {
			/* L40: */
			dky[i__ -1] = c__ * yh[i__ + jp1 * yh_dim1 -yh_offset] + s * dky[i__ -1];
		}
		/* L50: */
    }
    if (k == 0) {
		return 0;
    }
L55:
    i__1 = -(k);
    r__ = pow(h__, (double )i__1);
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L60: */
		dky[i__ -1] = r__ * dky[i__ -1];
    }
    return 0;
	
L80:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DINTDY-  K (=I1) illegal      \n");
#endif
    *iflag = -1;
    return 0;
L90:
#if (GMS_CUDA_DEBUG_ON) == 1
    fprintf(stderr, "DINTDY-  T (=R1) illegal      \n");

    fprintf(stderr, "      T not in interval TCUR - HU (= R1) to TCUR (=R2)      \n");
#endif
    *iflag = -2;
    return 0;
	/* ----------------------- END OF SUBROUTINE DINTDY ---------------------- */
} /* dintdy_ */

/* DECK DSTODA */
/* Subroutine */ 


template<typename Fex, typename Jex> 
__device__ int dstoda_(int * __restrict__ neq, 
                       double  * __restrict__ y, 
                       double  * __restrict__ yh, 
                       int * __restrict__ NOT_nyh, 
                       double  * __restrict__ yh1, 
                       double  * __restrict__ ewt, 
                       double  * __restrict__ savf, 
                       double  * __restrict__ acor, 
                       double  * __restrict__ wm, 
                       int * __restrict__ iwm, 
                       Fex f, 
                       Jex jac, 
                       struct cuLsodaCommonBlock * __restrict__ common)
{
    /* Initialized data */
	
	
    /* System generated locals */
    int yh_dim1 = 0;
	int yh_offset = 0;
	int i__1 = 0;
	int i__2 = 0;

    double  d__1 = 0.;
	double  d__2 = 0.;
	double  d__3 = 0.;
	
    /* Builtin functions */
	
    /* Local variables */
    int i__ = 0;
	int j = 0;
	int m = 0;
    double  r__ = 0.;
    int i1 = 0;
	int jb = 0;
    double  rh = 0.;
	double  rm = 0.;
	double  dm1 = 0.;
	double  dm2 = 0.;
    int lm1 = 0;
	int lm2 = 0;
    double  rh1 = 0.;
	double  rh2 = 0.;
	double  del = 0.;
	double  ddn = 0.;
    int ncf = 0;
    double  pdh = 0.;
	double  dsm = 0.;
	double  dup = 0.;
	double  exm1 = 0.;
	double  exm2 = 0.;
    int nqm1 = 0;
	int nqm2 = 0;
    double  dcon = 0.;
	double  delp = 0.;
    int lm1p1 = 0;
	int lm2p1 = 0;
    double  exdn = 0.;
	double  rhdn = 0.;
    int iret = 0;
    double  told = 0.;
	double  rhsm = 0.;
    int newq = 0;
    double  exsm = 0.;
	double  rhup = 0.;
	double  rate = 0.;
	double  exup = 0.;
	double  rh1it = 0.;
	double  alpha = 0.;
    int iredo = 0;
    double  pnorm = 0.;
	
    /* Parameter adjustments */
    //--neq;  //fixed
   // --y;
    yh_dim1 = *NOT_nyh;
    yh_offset = 1 + yh_dim1;
    //yh -= yh_offset;
    //--yh1;
    //--ewt;
    //--savf;
    //--acor;
    //--wm;
    //--iwm; 
	
    /* Function Body */
	/* ----------------------------------------------------------------------- */
	/* DSTODA performs one step of the integration of an initial value */
	/* problem for a system of ordinary differential equations. */
	/* Note: DSTODA is independent of the value of the iteration method */
	/* indicator MITER, when this is .ne. 0, and hence is independent */
	/* of the type of chord method used, or the Jacobian structure. */
	/* Communication with DSTODA is done with the following variables: */
	
	/* Y      = an array of length .ge. N used as the Y argument in */
	/*          all calls to F and JAC. */
	/* NEQ    = int array containing problem size in NEQ(1), and */
	/*          passed as the NEQ argument in all calls to F and JAC. */
	/* YH     = an NYH by LMAX array containing the dependent variables */
	/*          and their approximate scaled derivatives, where */
	/*          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate */
	/*          j-th derivative of y(i), scaled by H**j/factorial(j) */
	/*          (j = 0,1,...,NQ).  On entry for the first step, the first */
	/*          two columns of YH must be set from the initial values. */
	/* NYH    = a constant int .ge. N, the first dimension of YH. */
	/* YH1    = a one-dimensional array occupying the same space as YH. */
	/* EWT    = an array of length N containing multiplicative weights */
	/*          for local error measurements.  Local errors in y(i) are */
	/*          compared to 1.0/EWT(i) in various error tests. */
	/* SAVF   = an array of working storage, of length N. */
	/* ACOR   = a work array of length N, used for the accumulated */
	/*          corrections.  On a successful return, ACOR(i) contains */
	/*          the estimated one-step local error in y(i). */
	/* WM,IWM = real and int work arrays associated with matrix */
	/*          operations in chord iteration (MITER .ne. 0). */
	/* dprja_   = name of routine to evaluate and preprocess Jacobian matrix */
	/*          and P = I - H*EL0*Jac, if a chord method is being used. */
	/*          It also returns an estimate of norm(Jac) in PDNORM. */
	/* dsolsy_   = name of routine to solve linear system in chord iteration. */
	/* CCMAX  = maximum relative change in H*EL0 before dprja_ is called. */
	/* H      = the step size to be attempted on the next step. */
	/*          H is altered by the error control algorithm during the */
	/*          problem.  H can be either positive or negative, but its */
	/*          sign must remain constant throughout the problem. */
	/* HMIN   = the minimum absolute value of the step size H to be used. */
	/* HMXI   = inverse of the maximum absolute value of H to be used. */
	/*          HMXI = 0.0 is allowed and corresponds to an infinite HMAX. */
	/*          HMIN and HMXI may be changed at any time, but will not */
	/*          take effect until the next change of H is considered. */
	/* TN     = the independent variable. TN is updated on each step taken. */
	/* JSTART = an int used for input only, with the following */
	/*          values and meanings: */
	/*               0  perform the first step. */
	/*           .gt.0  take a new step continuing from the last. */
	/*              -1  take the next step with a new value of H, */
	/*                    N, METH, MITER, and/or matrix parameters. */
	/*              -2  take the next step with a new value of H, */
	/*                    but with other inputs unchanged. */
	/*          On return, JSTART is set to 1 to facilitate continuation. */
	/* KFLAG  = a completion code with the following meanings: */
	/*               0  the step was succesful. */
	/*              -1  the requested error could not be achieved. */
	/*              -2  corrector convergence could not be achieved. */
	/*              -3  fatal error in dprja_ or dsolsy_. */
	/*          A return with KFLAG = -1 or -2 means either */
	/*          ABS(H) = HMIN or 10 consecutive failures occurred. */
	/*          On a return with KFLAG negative, the values of TN and */
	/*          the YH array are as of the beginning of the last */
	/*          step, and H is the last step size attempted. */
	/* MAXORD = the maximum order of integration method to be allowed. */
	/* MAXCOR = the maximum number of corrector iterations allowed. */
	/* MSBP   = maximum number of steps between dprja_ calls (MITER .gt. 0). */
	/* MXNCF  = maximum number of convergence failures allowed. */
	/* METH   = current method. */
	/*          METH = 1 means Adams method (nonstiff) */
	/*          METH = 2 means BDF method (stiff) */
	/*          METH may be reset by DSTODA. */
	/* MITER  = corrector iteration method. */
	/*          MITER = 0 means functional iteration. */
	/*          MITER = JT .gt. 0 means a chord iteration corresponding */
	/*          to Jacobian type JT.  (The DLSODA/DLSODAR argument JT is */
	/*          communicated here as JTYP, but is not used in DSTODA */
	/*          except to load MITER following a method switch.) */
	/*          MITER may be reset by DSTODA. */
	/* N      = the number of first-order differential equations. */
	/* ----------------------------------------------------------------------- */
    kflag = 0;
    told = tn;
    ncf = 0;
    ierpj = 0;
    iersl = 0;
    jcur = 0;
    icf = 0;
    delp = 0.;
    if (jstart > 0) {
		goto L200;
    }
    if (jstart == -1) {
		goto L100;
    }
    if (jstart == -2) {
		goto L160;
    }
	/* ----------------------------------------------------------------------- */
	/* On the first call, the order is set to 1, and other variables are */
	/* initialized.  RMAX is the maximum ratio by which H can be increased */
	/* in a single step.  It is initially 1.E4 to compensate for the small */
	/* initial H, but then is normally equal to 10.  If a failure */
	/* occurs (in corrector convergence or error test), RMAX is set at 2 */
	/* for the next increase. */
	/* DCFODE is called to get the needed coefficients for both methods. */
	/* ----------------------------------------------------------------------- */
    lmax = maxord + 1;
    nq = 1;
    l = 2;
    ialth = 2;
    rmax = 1e4;
    rc = 0.;
    el0 = 1.;
    crate = .7;
    hold = h__;
    nslp = 0;
    ipup = miter;
    iret = 3;
	/* Initialize switching parameters.  METH = 1 is assumed initially. ----- */
    icount = 20;
    irflag = 0;
    pdest = 0.;
    pdlast = 0.;
    ratio = 5.;
    dcfode_(2, elco, tesco, common);
    for (i__ = 1; i__ <= 5; ++i__) {
		/* L10: */
		cm2[i__ - 1] = tesco[i__ * 3 - 2] * elco[i__ + 1 + i__ * 13 - 14];
    }
    dcfode_(1, elco, tesco, common);
    for (i__ = 1; i__ <= 12; ++i__) {
		/* L20: */
		cm1[i__ - 1] = tesco[i__ * 3 - 2] * elco[i__ + 1 + i__ * 13 - 14];
    }
    goto L150;
	/* ----------------------------------------------------------------------- */
	/* The following block handles preliminaries needed when JSTART = -1. */
	/* IPUP is set to MITER to force a matrix update. */
	/* If an order increase is about to be considered (IALTH = 1), */
	/* IALTH is reset to 2 to postpone consideration one more step. */
	/* If the caller has changed METH, DCFODE is called to reset */
	/* the coefficients of the method. */
	/* If H is to be changed, YH must be rescaled. */
	/* If H or METH is being changed, IALTH is reset to L = NQ + 1 */
	/* to prevent further changes in H for that many steps. */
	/* ----------------------------------------------------------------------- */
L100:
    ipup = miter;
    lmax = maxord + 1;
    if (ialth == 1) {
		ialth = 2;
    }
    if (meth == mused) {
		goto L160;
    }
    dcfode_(meth, elco, tesco, common);
    ialth = l;
    iret = 1;
	/* ----------------------------------------------------------------------- */
	/* The el vector and related constants are reset */
	/* whenever the order NQ is changed, or at the start of the problem. */
	/* ----------------------------------------------------------------------- */
L150:
    i__1 = l;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L155: */
		el[i__ - 1] = elco[i__ + nq * 13 - 14];
    }
    nqnyh = nq * *NOT_nyh;
    rc = rc * el[0] / el0;
    el0 = el[0];
    conit = .5 / (double )(nq + 2);
    switch (iret) {
		case 1:  goto L160;
		case 2:  goto L170;
		case 3:  goto L200;
    }
	/* ----------------------------------------------------------------------- */
	/* If H is being changed, the H ratio RH is checked against */
	/* RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to */
	/* L = NQ + 1 to prevent a change of H for that many steps, unless */
	/* forced by a convergence or error test failure. */
	/* ----------------------------------------------------------------------- */
L160:
    if (h__ == hold) {
		goto L200;
    }
    rh = h__ / hold;
    h__ = hold;
    iredo = 3;
    goto L175;
L170:
	/* Computing MAX */
    d__1 = rh, d__2 = hmin / fabs(h__);
    rh = max(d__1,d__2);
L175:
    rh = min(rh,rmax);
	/* Computing MAX */
    d__1 = 1., d__2 = fabs(h__) * hmxi * rh;
    rh /= max(d__1,d__2);
	/* ----------------------------------------------------------------------- */
	/* If METH = 1, also restrict the new step size by the stability region. */
	/* If this reduces H, set IRFLAG to 1 so that if there are roundoff */
	/* problems later, we can assume that is the cause of the trouble. */
	/* ----------------------------------------------------------------------- */
    if (meth == 2) {
		goto L178;
    }
    irflag = 0;
	/* Computing MAX */
    d__1 = fabs(h__) * pdlast;
    pdh = max(d__1,1e-6);
    if (rh * pdh * 1.00001 < sm1[nq - 1]) {
		goto L178;
    }
    rh = sm1[nq - 1] / pdh;
    irflag = 1;
L178:
    r__ = 1.;
    i__1 = l;
    for (j = 2; j <= i__1; ++j) {
		r__ *= rh;
		i__2 = n;
		for (i__ = 1; i__ <= i__2; ++i__) {
			/* L180: */
			yh[i__ + j * yh_dim1 -yh_offset] *= r__;
		}
    }
    h__ *= rh;
    rc *= rh;
    ialth = l;
    if (iredo == 0) {
		goto L690;
    }
	/* ----------------------------------------------------------------------- */
	/* This section computes the predicted values by effectively */
	/* multiplying the YH array by the Pascal triangle matrix. */
	/* RC is the ratio of new to old values of the coefficient  H*EL(1). */
	/* When RC differs from 1 by more than CCMAX, IPUP is set to MITER */
	/* to force dprja_ to be called, if a Jacobian is involved. */
	/* In any case, dprja_ is called at least every MSBP steps. */
	/* ----------------------------------------------------------------------- */
L200:
    if ((d__1 = rc - 1., fabs(d__1)) > ccmax) {
		ipup = miter;
    }
    if (nst >= nslp + msbp) {
		ipup = miter;
    }
    tn += h__;
    i1 = nqnyh + 1;
    i__2 = nq;
    for (jb = 1; jb <= i__2; ++jb) {
		i1 -= *NOT_nyh;
		/* DIR$ IVDEP */
		i__1 = nqnyh;
		for (i__ = i1; i__ <= i__1; ++i__) {
			/* L210: */
			yh1[i__ -1] += yh1[i__ + *NOT_nyh -1];
		}
		/* L215: */
    }
    pnorm = dmnorm_(&n, yh1, ewt, common);
	/* ----------------------------------------------------------------------- */
	/* Up to MAXCOR corrector iterations are taken.  A convergence test is */
	/* made on the RMS-norm of each correction, weighted by the error */
	/* weight vector EWT.  The sum of the corrections is accumulated in the */
	/* vector ACOR(i).  The YH array is not altered in the corrector loop. */
	/* ----------------------------------------------------------------------- */
L220:
    m = 0;
    rate = 0.;
    del = 0.;
    i__2 = n;
    for (i__ = 1; i__ <= i__2; ++i__) {
		/* L230: */
		y[i__ - 1] = yh[i__ + yh_dim1 -yh_offset]; //fixed y
    }
    f(neq, &tn, y, savf); //fixed neq y
    ++nfe;
    if (ipup <= 0) {
		goto L250;
    }
	/* ----------------------------------------------------------------------- */
	/* If indicated, the matrix P = I - H*EL(1)*J is reevaluated and */
	/* preprocessed before starting the corrector iteration.  IPUP is set */
	/* to 0 as an indicator that this has been done. */
	/* ----------------------------------------------------------------------- */
    dprja_(neq, y, &yh[yh_offset -yh_offset], NOT_nyh, ewt, acor, savf, wm, iwm, f, jac, common);  //fixed neq y
    ipup = 0;
    rc = 1.;
    nslp = nst;
    crate = .7;
    if (ierpj != 0) {
		goto L430;
    }
L250:
    i__2 = n;
    for (i__ = 1; i__ <= i__2; ++i__) {
		/* L260: */
		acor[i__ -1] = 0.;
    }
L270:
    if (miter != 0) {
		goto L350;
    }
	/* ----------------------------------------------------------------------- */
	/* In the case of functional iteration, update Y directly from */
	/* the result of the last function evaluation. */
	/* ----------------------------------------------------------------------- */
    i__2 = n;
    for (i__ = 1; i__ <= i__2; ++i__) {
		savf[i__ -1] = h__ * savf[i__ -1] - yh[i__ + (yh_dim1 << 1) -yh_offset];
		/* L290: */
		y[i__ - 1] = savf[i__ -1] - acor[i__ -1]; //fixed y
    }
    del = dmnorm_(&n, y, ewt, common);
    i__2 = n;
    for (i__ = 1; i__ <= i__2; ++i__) {
		y[i__ - 1] = yh[i__ + yh_dim1 -yh_offset] + el[0] * savf[i__ -1];
		/* L300: */
		acor[i__ -1] = savf[i__ -1];
    }
    goto L400;
	/* ----------------------------------------------------------------------- */
	/* In the case of the chord method, compute the corrector error, */
	/* and solve the linear system with that as right-hand side and */
	/* P as coefficient matrix. */
	/* ----------------------------------------------------------------------- */
L350:
    i__2 = n;
    for (i__ = 1; i__ <= i__2; ++i__) {
		/* L360: */
		y[i__ - 1] = h__ * savf[i__ -1] - (yh[i__ + (yh_dim1 << 1) -yh_offset] + acor[i__ -1]);
    }
    dsolsy_(wm, iwm, y, savf, common);
    if (iersl < 0) {
		goto L430;
    }
    if (iersl > 0) {
		goto L410;
    }
    del = dmnorm_(&n, y, ewt, common);
    i__2 = n;
    for (i__ = 1; i__ <= i__2; ++i__) {
		acor[i__ -1] += y[i__ - 1];
		/* L380: */
		y[i__ - 1] = yh[i__ + yh_dim1 -yh_offset] + el[0] * acor[i__ -1];
    }
	/* ----------------------------------------------------------------------- */
	/* Test for convergence.  If M .gt. 0, an estimate of the convergence */
	/* rate constant is stored in CRATE, and this is used in the test. */
	
	/* We first check for a change of iterates that is the size of */
	/* roundoff error.  If this occurs, the iteration has converged, and a */
	/* new rate estimate is not formed. */
	/* In all other cases, force at least two iterations to estimate a */
	/* local Lipschitz constant estimate for Adams methods. */
	/* On convergence, form PDEST = local maximum Lipschitz constant */
	/* estimate.  PDLAST is the most recent nonzero estimate. */
	/* ----------------------------------------------------------------------- */
L400:
    if (del <= pnorm * 100. * uround) {
		goto L450;
    }
    if (m == 0 && meth == 1) {
		goto L405;
    }
    if (m == 0) {
		goto L402;
    }
    rm = 1024.;
    if (del <= delp * 1024.) {
		rm = del / delp;
    }
    rate = max(rate,rm);
	/* Computing MAX */
    d__1 = crate * .2;
    crate = max(d__1,rm);
L402:
	/* Computing MIN */
    d__1 = 1., d__2 = crate * 1.5;
    dcon = del * min(d__1,d__2) / (tesco[nq * 3 - 2] * conit);
    if (dcon > 1.) {
		goto L405;
    }
	/* Computing MAX */
    d__2 = pdest, d__3 = rate / (d__1 = h__ * el[0]
								 , fabs(d__1));
    pdest = max(d__2,d__3);
    if (pdest != 0.) {
		pdlast = pdest;
    }
    goto L450;
L405:
    ++m;
    if (m == maxcor) {
		goto L410;
    }
    if (m >= 2 && del > delp * 2.) {
		goto L410;
    }
    delp = del;
    f(neq, &tn, y, savf); //fixed neq y
    ++nfe;
    goto L270;
	/* ----------------------------------------------------------------------- */
	/* The corrector iteration failed to converge. */
	/* If MITER .ne. 0 and the Jacobian is out of date, dprja_ is called for */
	/* the next try.  Otherwise the YH array is retracted to its values */
	/* before prediction, and H is reduced, if possible.  If H cannot be */
	/* reduced or MXNCF failures have occurred, exit with KFLAG = -2. */
	/* ----------------------------------------------------------------------- */
L410:
    if (miter == 0 || jcur == 1) {
		goto L430;
    }
    icf = 1;
    ipup = miter;
    goto L220;
L430:
    icf = 2;
    ++ncf;
    rmax = 2.;
    tn = told;
    i1 = nqnyh + 1;
    i__2 = nq;
    for (jb = 1; jb <= i__2; ++jb) {
		i1 -= *NOT_nyh;
		/* DIR$ IVDEP */
		i__1 = nqnyh;
		for (i__ = i1; i__ <= i__1; ++i__) {
			/* L440: */
			yh1[i__ -1] -= yh1[i__ + *NOT_nyh -1];
		}
		/* L445: */
    }
    if (ierpj < 0 || iersl < 0) {
		goto L680;
    }
    if (fabs(h__) <= hmin * 1.00001) {
		goto L670;
    }
    if (ncf == mxncf) {
		goto L670;
    }
    rh = .25;
    ipup = miter;
    iredo = 1;
    goto L170;
	/* ----------------------------------------------------------------------- */
	/* The corrector has converged.  JCUR is set to 0 */
	/* to signal that the Jacobian involved may need updating later. */
	/* The local error test is made and control passes to statement 500 */
	/* if it fails. */
	/* ----------------------------------------------------------------------- */
L450:
    jcur = 0;
    if (m == 0) {
		dsm = del / tesco[nq * 3 - 2];
    }
    if (m > 0) {
		dsm = dmnorm_(&n, acor, ewt, common) / tesco[nq * 3 - 2];
    }
    if (dsm > 1.) {
		goto L500;
    }
	/* ----------------------------------------------------------------------- */
	/* After a successful step, update the YH array. */
	/* Decrease ICOUNT by 1, and if it is -1, consider switching methods. */
	/* If a method switch is made, reset various parameters, */
	/* rescale the YH array, and exit.  If there is no switch, */
	/* consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1. */
	/* If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for */
	/* use in a possible order increase on the next step. */
	/* If a change in H is considered, an increase or decrease in order */
	/* by one is considered also.  A change in H is made only if it is by a */
	/* factor of at least 1.1.  If not, IALTH is set to 3 to prevent */
	/* testing for that many steps. */
	/* ----------------------------------------------------------------------- */
    kflag = 0;
    iredo = 0;
    ++nst;
    hu = h__;
    nqu = nq;
    mused = meth;
    i__2 = l;
    for (j = 1; j <= i__2; ++j) {
		i__1 = n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			/* L460: */
			yh[i__ + j * yh_dim1 -yh_offset] += el[j - 1] * acor[i__ -1];
		}
    }
    --icount;
    if (icount >= 0) {
		goto L488;
    }
    if (meth == 2) {
		goto L480;
    }
	/* ----------------------------------------------------------------------- */
	/* We are currently using an Adams method.  Consider switching to BDF. */
	/* If the current order is greater than 5, assume the problem is */
	/* not stiff, and skip this section. */
	/* If the Lipschitz constant and error estimate are not polluted */
	/* by roundoff, go to 470 and perform the usual test. */
	/* Otherwise, switch to the BDF methods if the last step was */
	/* restricted to insure stability (irflag = 1), and stay with Adams */
	/* method if not.  When switching to BDF with polluted error estimates, */
	/* in the absence of other information, double  the step size. */
	
	/* When the estimates are OK, we make the usual test by computing */
	/* the step size we could have (ideally) used on this step, */
	/* with the current (Adams) method, and also that for the BDF. */
	/* If NQ .gt. MXORDS, we consider changing to order MXORDS on switching. */
	/* Compare the two step sizes to decide whether to switch. */
	/* The step size advantage must be at least RATIO = 5 to switch. */
	/* ----------------------------------------------------------------------- */
    if (nq > 5) {
		goto L488;
    }
    if (dsm > pnorm * 100. * uround && pdest != 0.) {
		goto L470;
    }
    if (irflag == 0) {
		goto L488;
    }
    rh2 = 2.;
    nqm2 = min(nq,mxords);
    goto L478;
L470:
    exsm = 1. / (double )l;
    rh1 = 1. / (pow(dsm, exsm) * 1.2 + 1.2e-6);
    rh1it = rh1 * 2.;
    pdh = pdlast * fabs(h__);
    if (pdh * rh1 > 1e-5) {
		rh1it = sm1[nq - 1] / pdh;
    }
    rh1 = min(rh1,rh1it);
    if (nq <= mxords) {
		goto L474;
    }
    nqm2 = mxords;
    lm2 = mxords + 1;
    exm2 = 1. / (double )lm2;
    lm2p1 = lm2 + 1;
    dm2 = dmnorm_(&n, &yh[lm2p1 * yh_dim1 + 1 -yh_offset], ewt, common) / cm2[mxords - 1];
    rh2 = 1. / (pow(dm2, exm2) * 1.2 + 1.2e-6);
    goto L476;
L474:
    dm2 = dsm * (cm1[nq - 1] / cm2[nq - 1]
				 );
    rh2 = 1. / (pow(dm2, exsm) * 1.2 + 1.2e-6);
    nqm2 = nq;
L476:
    if (rh2 < ratio * rh1) {
		goto L488;
    }
	/* THE SWITCH TEST PASSED.  RESET RELEVANT QUANTITIES FOR BDF. ---------- */
L478:
    rh = rh2;
    icount = 20;
    meth = 2;
    miter = jtyp;
    pdlast = 0.;
    nq = nqm2;
    l = nq + 1;
    goto L170;
	/* ----------------------------------------------------------------------- */
	/* We are currently using a BDF method.  Consider switching to Adams. */
	/* Compute the step size we could have (ideally) used on this step, */
	/* with the current (BDF) method, and also that for the Adams. */
	/* If NQ .gt. MXORDN, we consider changing to order MXORDN on switching. */
	/* Compare the two step sizes to decide whether to switch. */
	/* The step size advantage must be at least 5/RATIO = 1 to switch. */
	/* If the step size for Adams would be so small as to cause */
	/* roundoff pollution, we stay with BDF. */
	/* ----------------------------------------------------------------------- */
L480:
    exsm = 1. / (double )l;
    if (mxordn >= nq) {
		goto L484;
    }
    nqm1 = mxordn;
    lm1 = mxordn + 1;
    exm1 = 1. / (double )lm1;
    lm1p1 = lm1 + 1;
    dm1 = dmnorm_(&n, &yh[lm1p1 * yh_dim1 + 1 -yh_offset], ewt, common) / cm1[mxordn - 1];
    rh1 = 1. / (pow(dm1, exm1) * 1.2 + 1.2e-6);
    goto L486;
L484:
    dm1 = dsm * (cm2[nq - 1] / cm1[nq - 1]
				 );
    rh1 = 1. / (pow(dm1, exsm) * 1.2 + 1.2e-6);
    nqm1 = nq;
    exm1 = exsm;
L486:
    rh1it = rh1 * 2.;
    pdh = pdnorm * fabs(h__);
    if (pdh * rh1 > 1e-5) {
		rh1it = sm1[nqm1 - 1] / pdh;
    }
    rh1 = min(rh1,rh1it);
    rh2 = 1. / (pow(dsm, exsm) * 1.2 + 1.2e-6);
    if (rh1 * ratio < rh2 * 5.) {
		goto L488;
    }
    alpha = max(.001,rh1);
    dm1 = pow(alpha, exm1) * dm1;
    if (dm1 <= uround * 1e3 * pnorm) {
		goto L488;
    }
	/* The switch test passed.  Reset relevant quantities for Adams. -------- */
    rh = rh1;
    icount = 20;
    meth = 1;
    miter = 0;
    pdlast = 0.;
    nq = nqm1;
    l = nq + 1;
    goto L170;
	
	/* No method switch is being made.  Do the usual step/order selection. -- */
L488:
    --ialth;
    if (ialth == 0) {
		goto L520;
    }
    if (ialth > 1) {
		goto L700;
    }
    if (l == lmax) {
		goto L700;
    }
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L490: */
		yh[i__ + lmax * yh_dim1 -yh_offset] = acor[i__ -1];
    }
    goto L700;
	/* ----------------------------------------------------------------------- */
	/* The error test failed.  KFLAG keeps track of multiple failures. */
	/* Restore TN and the YH array to their previous values, and prepare */
	/* to try the step again.  Compute the optimum step size for this or */
	/* one lower order.  After 2 or more failures, H is forced to decrease */
	/* by a factor of 0.2 or less. */
	/* ----------------------------------------------------------------------- */
L500:
    --kflag;
    tn = told;
    i1 = nqnyh + 1;
    i__1 = nq;
    for (jb = 1; jb <= i__1; ++jb) {
		i1 -= *NOT_nyh;
		/* DIR$ IVDEP */
		i__2 = nqnyh;
		for (i__ = i1; i__ <= i__2; ++i__) {
			/* L510: */
			yh1[i__ -1] -= yh1[i__ + *NOT_nyh -1];
		}
		/* L515: */
    }
    rmax = 2.;
    if (fabs(h__) <= hmin * 1.00001) {
		goto L660;
    }
    if (kflag <= -3) {
		goto L640;
    }
    iredo = 2;
    rhup = 0.;
    goto L540;
	/* ----------------------------------------------------------------------- */
	/* Regardless of the success or failure of the step, factors */
	/* RHDN, RHSM, and RHUP are computed, by which H could be multiplied */
	/* at order NQ - 1, order NQ, or order NQ + 1, respectively. */
	/* In the case of failure, RHUP = 0.0 to avoid an order increase. */
	/* The largest of these is determined and the new order chosen */
	/* accordingly.  If the order is to be increased, we compute one */
	/* additional scaled derivative. */
	/* ----------------------------------------------------------------------- */
L520:
    rhup = 0.;
    if (l == lmax) {
		goto L540;
    }
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L530: */
		savf[i__ -1] = acor[i__ -1] - yh[i__ + lmax * yh_dim1 -yh_offset];
    }
    dup = dmnorm_(&n, savf, ewt, common) / tesco[nq * 3 - 1];
    exup = (double )1. / (double )(l + 1);
    rhup = (double )1. / (pow(dup, exup) * (double )1.4 + (double )1.4e-6);
L540:
    exsm = (double )1. / l;
    rhsm = (double )1. / (pow(dsm, exsm) * (double )1.2 + (double )1.2e-6);
    rhdn = 0.;
    if (nq == 1) {
		goto L550;
    }
    ddn = dmnorm_(&n, &yh[l * yh_dim1 + 1 -yh_offset], ewt, common) / 
	tesco[nq * 3 - 3];
    exdn = (double )1. / (double )nq;
    rhdn = (double )1. / (pow(ddn, exdn) * (double )1.3 + (double )1.3e-6);
	/* If METH = 1, limit RH according to the stability region also. -------- */
L550:
    if (meth == 2) {
		goto L560;
    }
	/* Computing MAX */
    d__1 = fabs(h__) * pdlast;
    pdh = max(d__1,1e-6);
    if (l < lmax) {
		/* Computing MIN */
		d__1 = rhup, d__2 = sm1[l - 1] / pdh;
		rhup = min(d__1,d__2);
    }
	/* Computing MIN */
    d__1 = rhsm, d__2 = sm1[nq - 1] / pdh;
    rhsm = min(d__1,d__2);
    if (nq > 1) {
		/* Computing MIN */
		d__1 = rhdn, d__2 = sm1[nq - 2] / pdh;
		rhdn = min(d__1,d__2);
    }
    pdest = 0.;
L560:
    if (rhsm >= rhup) {
		goto L570;
    }
    if (rhup > rhdn) {
		goto L590;
    }
    goto L580;
L570:
    if (rhsm < rhdn) {
		goto L580;
    }
    newq = nq;
    rh = rhsm;
    goto L620;
L580:
    newq = nq - 1;
    rh = rhdn;
    if (kflag < 0 && rh > 1.) {
		rh = 1.;
    }
    goto L620;
L590:
    newq = l;
    rh = rhup;
    if (rh < 1.1) {
		goto L610;
    }
    r__ = el[l - 1] / (double )l;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L600: */
		yh[i__ + (newq + 1) * yh_dim1 -yh_offset] = acor[i__ -1] * r__;
    }
    goto L630;
L610:
    ialth = 3;
    goto L700;
	/* If METH = 1 and H is restricted by stability, bypass 10 percent test. */
L620:
    if (meth == 2) {
		goto L622;
    }
    if (rh * pdh * 1.00001 >= sm1[newq - 1]) {
		goto L625;
    }
L622:
    if (kflag == 0 && rh < 1.1) {
		goto L610;
    }
L625:
    if (kflag <= -2) {
		rh = min(rh,.2);
    }
	/* ----------------------------------------------------------------------- */
	/* If there is a change of order, reset NQ, L, and the coefficients. */
	/* In any case H is reset according to RH and the YH array is rescaled. */
	/* Then exit from 690 if the step was OK, or redo the step otherwise. */
	/* ----------------------------------------------------------------------- */
    if (newq == nq) {
		goto L170;
    }
L630:
    nq = newq;
    l = nq + 1;
    iret = 2;
    goto L150;
	/* ----------------------------------------------------------------------- */
	/* Control reaches this section if 3 or more failures have occured. */
	/* If 10 failures have occurred, exit with KFLAG = -1. */
	/* It is assumed that the derivatives that have accumulated in the */
	/* YH array have errors of the wrong order.  Hence the first */
	/* derivative is recomputed, and the order is set to 1.  Then */
	/* H is reduced by a factor of 10, and the step is retried, */
	/* until it succeeds or H reaches HMIN. */
	/* ----------------------------------------------------------------------- */
L640:
    if (kflag == -10) {
		goto L660;
    }
    rh = .1;
	/* Computing MAX */
    d__1 = hmin / fabs(h__);
    rh = max(d__1,rh);
    h__ *= rh;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L645: */
		y[i__ - 1] = yh[i__ + yh_dim1 -yh_offset];
    }
    f(neq, &tn, y, savf); //fixed neq
    ++nfe;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L650: */
		yh[i__ + (yh_dim1 << 1) -yh_offset] = h__ * savf[i__ -1];
    }
    ipup = miter;
    ialth = 5;
    if (nq == 1) {
		goto L200;
    }
    nq = 1;
    l = 2;
    iret = 3;
    goto L150;
	/* ----------------------------------------------------------------------- */
	/* All returns are made through this section.  H is saved in HOLD */
	/* to allow the caller to change H on the next step. */
	/* ----------------------------------------------------------------------- */
L660:
    kflag = -1;
    goto L720;
L670:
    kflag = -2;
    goto L720;
L680:
    kflag = -3;
    goto L720;
L690:
    rmax = 10.;
L700:
    r__ = 1. / tesco[nqu * 3 - 2];
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L710: */
		acor[i__ -1] *= r__;
    }
L720:
    hold = h__;
    jstart = 1;
    return 0;
	/* ----------------------- End of Subroutine DSTODA ---------------------- */
} /* dstoda_ */

/* DECK DCFODE */
/* Subroutine */ 
__device__ int dcfode_(int PARAM_meth, 
                       double  * __restrict__ DCFODE_elco, 
                       double  * __restrict__ DCFODE_tesco, 
                       struct cuLsodaCommonBlock *__restrict__ common)
{
    /* System generated locals */
    int i__1 = 0;
	
    /* Local variables */
     int i__ = 0;
	 int ib = 0;
     double  pc[12];
	 for (int bubb = 0; bubb < 12; bubb ++)
		{
			pc[bubb] = 0.;
		}
     int DCFODE_nq = 0;
     double  fnq = 0.;
     int nqm1 = 0;
	 int nqp1 = 0;
     double  ragq = 0.;
	 double  pint = 0.;
	 double  xpin = 0.;
	 double  fnqm1 = 0.;
	 double  agamq = 0.;
	 double  rqfac = 0.;
	 double  tsign = 0.;
	 double  rq1fac = 0.;
	
	/* ***BEGIN PROLOGUE  DCFODE */
	/* ***SUBSIDIARY */
	/* ***PURPOSE  Set ODE integrator coefficients. */
	/* ***TYPE      DOUBLE PRECISION (SCFODE-S, DCFODE-D) */
	/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
	/* ***DESCRIPTION */
	
	/*  DCFODE is called by the integrator routine to set coefficients */
	/*  needed there.  The coefficients for the current method, as */
	/*  given by the value of METH, are set for all orders and saved. */
	/*  The maximum order assumed here is 12 if METH = 1 and 5 if METH = 2. */
	/*  (A smaller value of the maximum order is also allowed.) */
	/*  DCFODE is called once at the beginning of the problem, */
	/*  and is not called again unless and until METH is changed. */
	
	/*  The ELCO array contains the basic method coefficients. */
	/*  The coefficients el(i), 1 .le. i .le. nq+1, for the method of */
	/*  order nq are stored in ELCO(i,nq).  They are given by a genetrating */
	/*  polynomial, i.e., */
	/*      l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq. */
	/*  For the implicit Adams methods, l(x) is given by */
	/*      dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) = 0. */
	/*  For the BDF methods, l(x) is given by */
	/*      l(x) = (x+1)*(x+2)* ... *(x+nq)/K, */
	/*  where         K = factorial(nq)*(1 + 1/2 + ... + 1/nq). */
	
	/*  The TESCO array contains test constants used for the */
	/*  local error test and the selection of step size and/or order. */
	/*  At order nq, TESCO(k,nq) is used for the selection of step */
	/*  size at order nq - 1 if k = 1, at order nq if k = 2, and at order */
	/*  nq + 1 if k = 3. */
	
	/* ***SEE ALSO  DLSODE */
	/* ***ROUTINES CALLED  (NONE) */
	/* ***REVISION HISTORY  (YYMMDD) */
	/*   791129  DATE WRITTEN */
	/*   890501  Modified prologue to SLATEC/LDOC format.  (FNF) */
	/*   890503  Minor cosmetic changes.  (FNF) */
	/*   930809  Renamed to allow single/double  precision versions. (ACH) */
	/* ***END PROLOGUE  DCFODE */
	/* **End */
	
	/* ***FIRST EXECUTABLE STATEMENT  DCFODE */
    /* Parameter adjustments */
    //DCFODE_tesco -= 4;
    //DCFODE_elco -= 14;
	
    /* Function Body */
    switch (PARAM_meth) {
		case 1:  goto L100;
		case 2:  goto L200;
    }
	
L100:
    DCFODE_elco[14 -14] = 1.;
    DCFODE_elco[15 -14] = 1.;
    DCFODE_tesco[4 -4] = 0.;
    DCFODE_tesco[5 -4] = 2.;
    DCFODE_tesco[7 -4] = 1.;
    DCFODE_tesco[39 -4] = 0.;
    pc[0] = (double )1.;
    rqfac =(double ) 1.;
    for (DCFODE_nq = 2; DCFODE_nq <= 12; ++DCFODE_nq) {
		/* ----------------------------------------------------------------------- */
		/* The PC array will contain the coefficients of the polynomial */
		/*     p(x) = (x+1)*(x+2)*...*(x+nq-1). */
		/* Initially, p(x) = 1. */
		/* ----------------------------------------------------------------------- */
		rq1fac = rqfac;
		rqfac /= (double )DCFODE_nq;
		nqm1 = DCFODE_nq - 1;
		fnqm1 = (double ) nqm1;
		nqp1 = DCFODE_nq + 1;
		/* Form coefficients of p(x)*(x+nq-1). ---------------------------------- */
		pc[DCFODE_nq - 1] = 0.;
		i__1 = nqm1;
		for (ib = 1; ib <= i__1; ++ib) {
			i__ = nqp1 - ib;
			/* L110: */
			pc[i__ - 1] = pc[i__ - 2] + fnqm1 * pc[i__ - 1];
		}
		pc[0] = fnqm1 * pc[0];
		/* Compute integral, -1 to 0, of p(x) and x*p(x). ----------------------- */
		pint = pc[0];
		xpin = pc[0] / 2.;
		tsign = 1.;
		i__1 = DCFODE_nq;
		for (i__ = 2; i__ <= i__1; ++i__) {
			tsign = -tsign;
			pint += tsign * pc[i__ - 1] / (double )i__;
			/* L120: */
			xpin += tsign * pc[i__ - 1] / (double )(i__ + 1);
		}
		/* Store coefficients in ELCO and TESCO. -------------------------------- */
		DCFODE_elco[DCFODE_nq * 13 + 1 -14] = pint * rq1fac;
		DCFODE_elco[DCFODE_nq * 13 + 2 -14] = 1.;
		i__1 = DCFODE_nq;
		for (i__ = 2; i__ <= i__1; ++i__) {
			/* L130: */
			DCFODE_elco[i__ + 1 + DCFODE_nq * 13 -14] = rq1fac * pc[i__ - 1] / (double )i__;
		}
		agamq = rqfac * xpin;
		ragq = 1. / agamq;
		DCFODE_tesco[DCFODE_nq * 3 + 2 -4] = ragq;
		if (DCFODE_nq < 12) {
			DCFODE_tesco[nqp1 * 3 + 1 -4] = ragq * rqfac / (double )nqp1;
		}
		DCFODE_tesco[nqm1 * 3 + 3 -4] = ragq;
		/* L140: */
    }
    return 0;
	
L200:
    pc[0] = 1.;
    rq1fac = 1.;
    for (DCFODE_nq = 1; DCFODE_nq <= 5; ++DCFODE_nq) {
		/* ----------------------------------------------------------------------- */
		/* The PC array will contain the coefficients of the polynomial */
		/*     p(x) = (x+1)*(x+2)*...*(x+nq). */
		/* Initially, p(x) = 1. */
		/* ----------------------------------------------------------------------- */
		fnq = (double ) DCFODE_nq;
		nqp1 = DCFODE_nq + 1;
		/* Form coefficients of p(x)*(x+nq). ------------------------------------ */
		pc[nqp1 - 1] = 0.;
		i__1 = DCFODE_nq;
		for (ib = 1; ib <= i__1; ++ib) {
			i__ = DCFODE_nq + 2 - ib;
			/* L210: */
			pc[i__ - 1] = pc[i__ - 2] + fnq * pc[i__ - 1];
		}
		pc[0] = fnq * pc[0];
		/* Store coefficients in ELCO and TESCO. -------------------------------- */
		i__1 = nqp1;
		for (i__ = 1; i__ <= i__1; ++i__) {
			/* L220: */
			DCFODE_elco[i__ + DCFODE_nq * 13 -14] = pc[i__ - 1] / pc[1];
		}
		DCFODE_elco[DCFODE_nq * 13 + 2 -14] = 1.;
		DCFODE_tesco[DCFODE_nq * 3 + 1 -4] = rq1fac;
		DCFODE_tesco[DCFODE_nq * 3 + 2 -4] = ((double )nqp1) / DCFODE_elco[DCFODE_nq * 13 + 1 -14];
		DCFODE_tesco[DCFODE_nq * 3 + 3 -4] = ((double )(DCFODE_nq + 2)) / DCFODE_elco[DCFODE_nq * 13 + 1 -14];
		rq1fac /= fnq;
		/* L230: */
    }
    return 0;
	/* ----------------------- END OF SUBROUTINE DCFODE ---------------------- */
} /* dcfode_ */

/* DECK DPRJA */
/* Subroutine */
#ifdef use_export
export
#endif 
template<typename Fex, typename Jex> 
__device__ int dprja_(int *  __restrict__ neq, 
                      double  * __restrict__ y,
                      double  * __restrict__ yh, 
                      int * __restrict__ NOT_nyh,
                      double  * __restrict__ ewt, 
                      double  * __restrict__ ftem, 
                      double  * __restrict__ savf, 
                      double  * __restrict__ wm, 
                      int * __restrict__ iwm, 
                      Fex f, 
                      Jex jac, 
                      struct cuLsodaCommonBlock * __restrict__ common)
{
    /* System generated locals */
    int yh_dim1, yh_offset, i__1, i__2, i__3, i__4;
    double  d__1 = 0.;
	double  d__2 = 0.;
	
    /* Local variables */
    int i__ = 0;
	int j;
    double  r__;
    int i1, i2, j1;
    double  r0;
    int ii = 0;
	int jj = 0;
	int ml = 0;
	int mu = 0;
    double  yi = 0.;
	double  yj = 0.;
	double  hl0;
    int ml3 = 0;
	int np1 = 0;
    double  fac;
    int mba = 0;
	int ier = 0;
    double  con = 0.;
	double  yjj;
    int meb1 = 0;
	int lenp = 0;
    double  srur;
	int mband = 0;
	int meband = 0;
	
	/* ----------------------------------------------------------------------- */
	/* DPRJA is called by DSTODA to compute and process the matrix */
	/* P = I - H*EL(1)*J , where J is an approximation to the Jacobian. */
	/* Here J is computed by the user-supplied routine JAC if */
	/* MITER = 1 or 4 or by finite differencing if MITER = 2 or 5. */
	/* J, scaled by -H*EL(1), is stored in WM.  Then the norm of J (the */
	/* matrix norm consistent with the weighted max-norm on vectors given */
	/* by DMNORM) is computed, and J is overwritten by P.  P is then */
	/* subjected to LU decomposition in preparation for later solution */
	/* of linear systems with P as coefficient matrix.  This is done */
	/* by DGEFA if MITER = 1 or 2, and by DGBFA if MITER = 4 or 5. */
	
	/* In addition to variables described previously, communication */
	/* with DPRJA uses the following: */
	/* Y     = array containing predicted values on entry. */
	/* FTEM  = work array of length N (ACOR in DSTODA). */
	/* SAVF  = array containing f evaluated at predicted y. */
	/* WM    = real work space for matrices.  On output it contains the */
	/*         LU decomposition of P. */
	/*         Storage of matrix elements starts at WM(3). */
	/*         WM also contains the following matrix-related data: */
	/*         WM(1) = SQRT(UROUND), used in numerical Jacobian increments. */
	/* IWM   = int work space containing pivot information, starting at */
	/*         IWM(21).   IWM also contains the band parameters */
	/*         ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5. */
	/* EL0   = EL(1) (input). */
	/* PDNORM= norm of Jacobian matrix. (Output). */
	/* IERPJ = output error flag,  = 0 if no trouble, .gt. 0 if */
	/*         P matrix found to be singular. */
	/* JCUR  = output flag = 1 to indicate that the Jacobian matrix */
	/*         (or approximation) is now current. */
	/* This routine also uses the Common variables EL0, H, TN, UROUND, */
	/* MITER, N, NFE, and NJE. */
	/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
   // --neq;
    //--y;
    yh_dim1 = *NOT_nyh;
    yh_offset = 1 + yh_dim1;
    //yh -= yh_offset;
    //--ewt;
    //--ftem;
    //--savf;
    //--wm;
    //--iwm;
	
    /* Function Body */
    ++nje;
    ierpj = 0;
    jcur = 1;
    hl0 = h__ * el0;
    switch (miter) {
		case 1:  goto L100;
		case 2:  goto L200;
		case 3:  goto L300;
		case 4:  goto L400;
		case 5:  goto L500;
    }
	/* If MITER = 1, call JAC and multiply by scalar. ----------------------- */
L100:
    lenp = n * n;
    i__1 = lenp;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L110: */
		wm[i__ + 2 -1] = 0.;
    }
    jac(neq, &tn, y, 0, 0, &wm[3 -1], n); //fixed neq
    con = -hl0;
    i__1 = lenp;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L120: */
		wm[i__ + 2 -1] *= con;
    }
    goto L240;
	/* If MITER = 2, make N calls to F to approximate J. -------------------- */
L200:
    fac = dmnorm_(&n, savf, ewt, common);
    r0 = fabs(h__) * 1e3 * uround * ((double )n) * fac;
    if (r0 == 0.) {
		r0 = 1.;
    }
    srur = wm[0];
    j1 = 2;
    i__1 = n;
    for (j = 1; j <= i__1; ++j) {
		yj = y[j - 1];
		/* Computing MAX */
		d__1 = srur * fabs(yj), d__2 = r0 / ewt[j -1];
		r__ = max(d__1,d__2);
		y[j - 1] += r__;
		fac = -hl0 / r__;
		f(neq, &tn, y, ftem); //fixed neq
		i__2 = n;
		for (i__ = 1; i__ <= i__2; ++i__) {
			/* L220: */
			wm[i__ + j1 -1] = (ftem[i__ -1] - savf[i__ -1]) * fac;
		}
		y[j -1] = yj;
		j1 += n;
		/* L230: */
    }
    nfe += n;
L240:
	/* Compute norm of Jacobian. -------------------------------------------- */
    pdnorm = dfnorm_(&n, &wm[3 -1], ewt, common) / fabs(hl0);
	/* Add identity matrix. ------------------------------------------------- */
    j = 3;
    np1 = n + 1;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		wm[j -1] += 1.;
		/* L250: */
		j += np1;
    }
	/* Do LU decomposition on P. -------------------------------------------- */
    dgefa_(&wm[3 -1], &n, &n, &iwm[21 -1], &ier, common);
    if (ier != 0) {
		ierpj = 1;
    }
    return 0;
	/* Dummy block only, since MITER is never 3 in this routine. ------------ */
L300:
    return 0;
	/* If MITER = 4, call JAC and multiply by scalar. ----------------------- */
L400:
    ml = iwm[0];
    mu = iwm[1];
    ml3 = ml + 3;
    mband = ml + mu + 1;
    meband = mband + ml;
    lenp = meband * n;
    i__1 = lenp;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L410: */
		wm[i__ + 2 -1] = 0.;
    }
    jac(neq, &tn, y, ml, mu, &wm[ml3 -1], meband); //fixed neq
    con = -hl0;
    i__1 = lenp;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L420: */
		wm[i__ + 2 -1] *= con;
    }
    goto L570;
	/* If MITER = 5, make MBAND calls to F to approximate J. ---------------- */
L500:
    ml = iwm[0];
    mu = iwm[1];
    mband = ml + mu + 1;
    mba = min(mband,n);
    meband = mband + ml;
    meb1 = meband - 1;
    srur = wm[0];
    fac = dmnorm_(&n, savf, ewt, common);
    r0 = fabs(h__) * 1e3 * uround * n * fac;
    if (r0 == 0.) {
		r0 = 1.;
    }
    i__1 = mba;
    for (j = 1; j <= i__1; ++j) {
		i__2 = n;
		i__3 = mband;
		for (i__ = j; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3) {
			yi = y[i__ -1];
			/* Computing MAX */
			d__1 = srur * fabs(yi), d__2 = r0 / ewt[i__ -1];
			r__ = max(d__1,d__2);
			/* L530: */
			y[i__ - 1] += r__;
		}
		f(neq, &tn, y, ftem); //fixed neq
		i__3 = n;
		i__2 = mband;
		for (jj = j; i__2 < 0 ? jj >= i__3 : jj <= i__3; jj += i__2) {
			y[jj - 1] = yh[jj + yh_dim1 -yh_offset];
			yjj = y[jj - 1];
			/* Computing MAX */
			d__1 = srur * fabs(yjj), d__2 = r0 / ewt[jj -1];
			r__ = max(d__1,d__2);
			fac = -hl0 / r__;
			/* Computing MAX */
			i__4 = jj - mu;
			i1 = max(i__4,1);
			/* Computing MIN */
			i__4 = jj + ml;
			i2 = min(i__4,n);
			ii = jj * meb1 - ml + 2;
			i__4 = i2;
			for (i__ = i1; i__ <= i__4; ++i__) {
				/* L540: */
				wm[ii + i__ -1] = (ftem[i__ -1] - savf[i__ -1]) * fac;
			}
			/* L550: */
		}
		/* L560: */
    }
    nfe += mba;
L570:
	/* Compute norm of Jacobian. -------------------------------------------- */
    pdnorm = dbnorm_(&n, &wm[ml + 3 -1], &meband, &ml, &mu, ewt, common) / fabs(hl0);
	/* Add identity matrix. ------------------------------------------------- */
    ii = mband + 2;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		wm[ii -1] += 1.;
		/* L580: */
		ii += meband;
    }
	/* Do LU decomposition of P. -------------------------------------------- */
    dgbfa_(&wm[3 -1], &meband, &n, &ml, &mu, &iwm[21 -1], &ier, common);
    if (ier != 0) {
		ierpj = 1;
    }
    return 0;
	/* ----------------------- End of Subroutine DPRJA ----------------------- */
} /* dprja_ */

/* DECK DSOLSY */
/* Subroutine */ 
__device__ int dsolsy_(double  * __restrict__ wm, 
                       int  * __restrict__ iwm, 
                       double  * __restrict__ x, 
                       double  * __restrict__ tem, 
                       struct cuLsodaCommonBlock * __restrict__ common)
{
    /* System generated locals */
    int i__1 = 0;
	
    /* Local variables */
    int i__ = 0;
    double  r__ = 0.;
	double  di = 0.;
    int ml = 0;
	int mu = 0;
    double  hl0 = 0.;
	double  phl0 = 0.;
    int meband = 0;
	
	/* ***BEGIN PROLOGUE  DSOLSY */
	/* ***SUBSIDIARY */
	/* ***PURPOSE  ODEPACK linear system solver. */
	/* ***TYPE      DOUBLE PRECISION (SSOLSY-S, DSOLSY-D) */
	/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
	/* ***DESCRIPTION */
	
	/*  This routine manages the solution of the linear system arising from */
	/*  a chord iteration.  It is called if MITER .ne. 0. */
	/*  If MITER is 1 or 2, it calls DGESL to accomplish this. */
	/*  If MITER = 3 it updates the coefficient h*EL0 in the diagonal */
	/*  matrix, and then computes the solution. */
	/*  If MITER is 4 or 5, it calls DGBSL. */
	/*  Communication with DSOLSY uses the following variables: */
	/*  WM    = real work space containing the inverse diagonal matrix if */
	/*          MITER = 3 and the LU decomposition of the matrix otherwise. */
	/*          Storage of matrix elements starts at WM(3). */
	/*          WM also contains the following matrix-related data: */
	/*          WM(1) = SQRT(UROUND) (not used here), */
	/*          WM(2) = HL0, the previous value of h*EL0, used if MITER = 3. */
	/*  IWM   = int work space containing pivot information, starting at */
	/*          IWM(21), if MITER is 1, 2, 4, or 5.  IWM also contains band */
	/*          parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5. */
	/*  X     = the right-hand side vector on input, and the solution vector */
	/*          on output, of length N. */
	/*  TEM   = vector of work space of length N, not used in this version. */
	/*  IERSL = output flag (in COMMON).  IERSL = 0 if no trouble occurred. */
	/*          IERSL = 1 if a singular matrix arose with MITER = 3. */
	/*  This routine also uses the COMMON variables EL0, H, MITER, and N. */
	
	/* ***SEE ALSO  DLSODE */
	/* ***ROUTINES CALLED  DGBSL, DGESL */
	/* ***COMMON BLOCKS    DLS001 */
	/* ***REVISION HISTORY  (YYMMDD) */
	/*   791129  DATE WRITTEN */
	/*   890501  Modified prologue to SLATEC/LDOC format.  (FNF) */
	/*   890503  Minor cosmetic changes.  (FNF) */
	/*   930809  Renamed to allow single/double  precision versions. (ACH) */
	/*   010418  Reduced size of Common block /DLS001/. (ACH) */
	/*   031105  Restored 'own' variables to Common block /DLS001/, to */
	/*           enable interrupt/restart feature. (ACH) */
	/* ***END PROLOGUE  DSOLSY */
	/* **End */
	
	/* ***FIRST EXECUTABLE STATEMENT  DSOLSY */
    /* Parameter adjustments */
    //--tem;
    //--x;
    //--iwm;
    //--wm;
	
    /* Function Body */
    iersl = 0;
    switch (miter) {
		case 1:  goto L100;
		case 2:  goto L100;
		case 3:  goto L300;
		case 4:  goto L400;
		case 5:  goto L400;
    }
L100:
    dgesl_(&wm[3 -1], &n, &n, &iwm[21 -1], x, 0, common);
    return 0;
	
L300:
    phl0 = wm[1];
    hl0 = h__ * el0;
    wm[1] = hl0;
    if (hl0 == phl0) {
		goto L330;
    }
    r__ = hl0 / phl0;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		di = 1. - r__ * (1. - 1. / wm[i__ + 2 -1]);
		if (fabs(di) == 0.) {
			goto L390;
		}
		/* L320: */
		wm[i__ + 2 -1] = 1. / di;
    }
L330:
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L340: */
		x[i__ -1] = wm[i__ + 2 -1] * x[i__ -1];
    }
    return 0;
L390:
    iersl = 1;
    return 0;
	
L400:
    ml = iwm[0];
    mu = iwm[1];
    meband = (ml << 1) + mu + 1;
    dgbsl_(&wm[3 -1], &meband, &n, &ml, &mu, &iwm[21 -1], x, 0, common);
    return 0;
	/* ----------------------- END OF SUBROUTINE DSOLSY ---------------------- */
} /* dsolsy_ */

/* DECK DEWSET */
/* Subroutine */ 
__device__ int dewset_(int * __restrict__ PARAM_n, 
                       int * __restrict__ itol, 
                       double  * __restrict__ rtol,  
                       double  * __restrict__ atol, 
                       double  * __restrict__ ycur, 
                       double  * __restrict__ ewt, 
                       struct cuLsodaCommonBlock * __restrict__ common)
{
    /* System generated locals */
    int i__1 = 0;
    double  d__1 = 0.;
	
    /* Local variables */
	int i__ = 0;
	
	/* ***BEGIN PROLOGUE  DEWSET */
	/* ***SUBSIDIARY */
	/* ***PURPOSE  Set error weight vector. */
	/* ***TYPE      DOUBLE PRECISION (SEWSET-S, DEWSET-D) */
	/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
	/* ***DESCRIPTION */
	
	/*  This subroutine sets the error weight vector EWT according to */
	/*      EWT(i) = RTOL(i)*ABS(YCUR(i)) + ATOL(i),  i = 1,...,N, */
	/*  with the subscript on RTOL and/or ATOL possibly replaced by 1 above, */
	/*  depending on the value of ITOL. */
	
	/* ***SEE ALSO  DLSODE */
	/* ***ROUTINES CALLED  (NONE) */
	/* ***REVISION HISTORY  (YYMMDD) */
	/*   791129  DATE WRITTEN */
	/*   890501  Modified prologue to SLATEC/LDOC format.  (FNF) */
	/*   890503  Minor cosmetic changes.  (FNF) */
	/*   930809  Renamed to allow single/double  precision versions. (ACH) */
	/* ***END PROLOGUE  DEWSET */
	/* **End */
	
	/* ***FIRST EXECUTABLE STATEMENT  DEWSET */
    /* Parameter adjustments */
    //--ewt;
    //--ycur;
    //--rtol;
    //--atol;
	
    /* Function Body */
    switch (*itol) {
		case 1:  goto L10;
		case 2:  goto L20;
		case 3:  goto L30;
		case 4:  goto L40;
    }
L10:
    i__1 = *PARAM_n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L15: */
		ewt[i__ -1] = rtol[0] * (d__1 = ycur[i__ -1], fabs(d__1)) + atol[0];
    }
    return 0;
L20:
    i__1 = *PARAM_n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L25: */
		ewt[i__ -1] = rtol[0] * (d__1 = ycur[i__ -1], fabs(d__1)) + atol[i__ -1];
    }
    return 0;
L30:
    i__1 = *PARAM_n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L35: */
		ewt[i__ -1] = rtol[i__ - 1] * (d__1 = ycur[i__ -1], fabs(d__1)) + atol[0];
    }
    return 0;
L40:
    i__1 = *PARAM_n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L45: */
		ewt[i__ -1] = rtol[i__ - 1] * (d__1 = ycur[i__ -1], fabs(d__1)) + atol[i__ -1];
    }
    return 0;
	/* ----------------------- END OF SUBROUTINE DEWSET ---------------------- */
} /* dewset_ */

/* DECK DMNORM */
__device__ double  dmnorm_(int * __restrict__ PARAM_n, double  * __restrict__ v, double  * __restrict__ w, struct cuLsodaCommonBlock * __restrict__ common)
{
    /* System generated locals */
    int i__1 = 0;
    double  ret_val = 0.;
	double  d__1 = 0.;
	double  d__2 = 0.;
	double  d__3 = 0.;
	
    /* Local variables */
     int i__ = 0;
     double  vm = 0.;
	
	/* ----------------------------------------------------------------------- */
	/* This function routine computes the weighted max-norm */
	/* of the vector of length N contained in the array V, with weights */
	/* contained in the array w of length N: */
	/*   DMNORM = MAX(i=1,...,N) ABS(V(i))*W(i) */
	/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    //--w;
    //--v;
	
    /* Function Body */
    vm = 0.;
    i__1 = *PARAM_n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L10: */
		/* Computing MAX */
		d__2 = vm, d__3 = (d__1 = v[i__ -1], fabs(d__1)) * w[i__ -1];
		vm = max(d__2,d__3);
    }
    ret_val = vm;
    return ret_val;
	/* ----------------------- End of Function DMNORM ------------------------ */
} /* dmnorm_ */

/* DECK DFNORM */
__device__ double  dfnorm_(int * __restrict__ PARAM_n, double  * __restrict__ a, double  * __restrict__ w, struct cuLsodaCommonBlock * __restrict__ common)
{
    /* System generated locals */
    int a_dim1 = 0;
	int a_offset = 0;
	int i__1 = 0;
	int i__2 = 0;
    double  ret_val = 0.;
	double  d__1 = 0.;
	double  d__2 = 0.;
	
    /* Local variables */
     int i__ = 0;
	 int j = 0;
     double  an = 0.;
	 double  sum = 0.;
	
	/* ----------------------------------------------------------------------- */
	/* This function computes the norm of a full N by N matrix, */
	/* stored in the array A, that is consistent with the weighted max-norm */
	/* on vectors, with weights stored in the array W: */
	/*   DFNORM = MAX(i=1,...,N) ( W(i) * Sum(j=1,...,N) ABS(a(i,j))/W(j) ) */
	/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    //--w;
    a_dim1 = *PARAM_n;
    a_offset = 1 + a_dim1;
    //a -= a_offset;
	
    /* Function Body */
    an = 0.;
    i__1 = *PARAM_n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		sum = 0.;
		i__2 = *PARAM_n;
		for (j = 1; j <= i__2; ++j) {
			/* L10: */
			sum += (d__1 = a[i__ + j * a_dim1 -a_offset], fabs(d__1)) / w[j -1];
		}
		/* Computing MAX */
		d__1 = an, d__2 = sum * w[i__ -1];
		an = max(d__1,d__2);
		/* L20: */
    }
    ret_val = an;
    return ret_val;
	/* ----------------------- End of Function DFNORM ------------------------ */
} /* dfnorm_ */

/* DECK DBNORM */
__device__ double  dbnorm_(int * __restrict__ PARAM_n, double  * __restrict__ a, int * __restrict__ nra, 
                        int * __restrict__ ml, int * __restrict__ mu, double  * __restrict__ w, struct cuLsodaCommonBlock * __restrict__ common)
{
 
	  /* System generated locals */
    int a_dim1 = 0;
	int a_offset = 0;
	int i__1 = 0;
	int i__2 = 0;
    double  ret_val = 0.;
	double  d__1 = 0.;
	double  d__2 = 0.;
	
    /* Local variables */
     int i__ = 0;
	 int j = 0;
     double  an = 0.;
	 double  sum = 0.;
     int i1 = 0;
     int jhi = 0;
	 int jlo = 0;

	
	/* ----------------------------------------------------------------------- */
	/* This function computes the norm of a banded N by N matrix, */
	/* stored in the array A, that is consistent with the weighted max-norm */
	/* on vectors, with weights stored in the array W. */
	/* ML and MU are the lower and upper half-bandwidths of the matrix. */
	/* NRA is the first dimension of the A array, NRA .ge. ML+MU+1. */
	/* In terms of the matrix elements a(i,j), the norm is given by: */
	/*   DBNORM = MAX(i=1,...,N) ( W(i) * Sum(j=1,...,N) ABS(a(i,j))/W(j) ) */
	/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    //--w;
    a_dim1 = *nra;
    a_offset = 1 + a_dim1;
    //a -= a_offset;
	
    /* Function Body */
    an = 0.;
    i__1 = *PARAM_n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		sum = 0.;
		i1 = i__ + *mu + 1;
		/* Computing MAX */
		i__2 = i__ - *ml;
		jlo = max(i__2,1);
		/* Computing MIN */
		i__2 = i__ + *mu;
		jhi = min(i__2,*PARAM_n);
		i__2 = jhi;
		for (j = jlo; j <= i__2; ++j) {
			/* L10: */
			sum += (d__1 = a[i1 - j + j * a_dim1 -a_offset], fabs(d__1)) / w[j -1];
		}
		/* Computing MAX */
		d__1 = an, d__2 = sum * w[i__ -1];
		an = max(d__1,d__2);
		/* L20: */
    }
    ret_val = an;
    return ret_val;
	/* ----------------------- End of Function DBNORM ------------------------ */
} /* dbnorm_ */

/* DECK DSRCMA */
/* Subroutine */ 
__device__ int dsrcma_(double  * __restrict__ rsav, int * __restrict__ isav, int * __restrict__job, struct cuLsodaCommonBlock * __restrict__ common)
{
    /* Initialized data */
	
     int lenrls = 218;
     int lenils = 37;
     int lenrla = 22;
     int lenila = 9;
	
    /* System generated locals */
    int i__1 = 0;
	
    /* Local variables */
     int i__ = 0;
	
	/* ----------------------------------------------------------------------- */
	/* This routine saves or restores (depending on JOB) the contents of */
	/* the Common blocks DLS001, DLSA01, which are used */
	/* internally by one or more ODEPACK solvers. */
	
	/* RSAV = real array of length 240 or more. */
	/* ISAV = int array of length 46 or more. */
	/* JOB  = flag indicating to save or restore the Common blocks: */
	/*        JOB  = 1 if Common is to be saved (written to RSAV/ISAV) */
	/*        JOB  = 2 if Common is to be restored (read from RSAV/ISAV) */
	/*        A call with JOB = 2 presumes a prior call with JOB = 1. */
	/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    //--isav;
    //--rsav;
	
    /* Function Body */
	
    if (*job == 2) {
		goto L100;
    }
    i__1 = lenrls;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L10: */
		rsav[i__ -1] = rls[i__ - 1];
    }
    i__1 = lenrla;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L15: */
		rsav[lenrls + i__ -1] = rlsa[i__ - 1];
    }
	
    i__1 = lenils;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L20: */
		isav[i__ -1] = ils[i__ - 1];
    }
    i__1 = lenila;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L25: */
		isav[lenils + i__ -1] = ilsa[i__ - 1];
    }
	
    return 0;
	
L100:
    i__1 = lenrls;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L110: */
		rls[i__ - 1] = rsav[i__ -1];
    }
    i__1 = lenrla;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L115: */
		rlsa[i__ - 1] = rsav[lenrls + i__ -1];
    }
	
    i__1 = lenils;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L120: */
		ils[i__ - 1] = isav[i__ -1];
    }
    i__1 = lenila;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* L125: */
		ilsa[i__ - 1] = isav[lenils + i__ -1];
    }
	
    return 0;
	/* ----------------------- End of Subroutine DSRCMA ---------------------- */
} /* dsrcma_ */

/* DECK DGEFA */
/* Subroutine */ 
__device__ int dgefa_(double  * __restrict__ a, int * __restrict__ lda, int * __restrict__ PARAM_n, 
                     int * __restrict__ ipvt, int * __restrict__ info, struct cuLsodaCommonBlock * __restrict__ common)
{
    /* System generated locals */
    int a_dim1 = 0;
	int a_offset = 0;
	int i__1 = 0;
	int i__2 = 0;
	int i__3 = 0;
	
    /* Local variables */
     int j = 0;
	 int k = 0;
	 int DGEFA_l = 0;
     double  t = 0.;
     int kp1 = 0;
	 int nm1 = 0;
	
	
	/* ***BEGIN PROLOGUE  DGEFA */
	/* ***PURPOSE  Factor a matrix using Gaussian elimination. */
	/* ***CATEGORY  D2A1 */
	/* ***TYPE      DOUBLE PRECISION (SGEFA-S, DGEFA-D, CGEFA-C) */
	/* ***KEYWORDS  GENERAL MATRIX, LINEAR ALGEBRA, LINPACK, */
	/*             MATRIX FACTORIZATION */
	/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
	/* ***DESCRIPTION */
	
	/*     DGEFA factors a double  precision matrix by Gaussian elimination. */
	
	/*     DGEFA is usually called by DGECO, but it can be called */
	/*     directly with a saving in time if  RCOND  is not needed. */
	/*     (Time for DGECO) = (1 + 9/N)*(Time for DGEFA) . */
	
	/*     On Entry */
	
	/*        A       DOUBLE PRECISION(LDA, N) */
	/*                the matrix to be factored. */
	
	/*        LDA     int */
	/*                the leading dimension of the array  A . */
	
	/*        N       int */
	/*                the order of the matrix  A . */
	
	/*     On Return */
	
	/*        A       an upper triangular matrix and the multipliers */
	/*                which were used to obtain it. */
	/*                The factorization can be written  A = L*U  where */
	/*                L  is a product of permutation and unit lower */
	/*                triangular matrices and  U  is upper triangular. */
	
	/*        IPVT    int(N) */
	/*                an int vector of pivot indices. */
	
	/*        INFO    int */
	/*                = 0  normal value. */
	/*                = K  if  U(K,K) .EQ. 0.0 .  This is not an error */
	/*                     condition for this subroutine, but it does */
	/*                     indicate that DGESL or DGEDI will divide by zero */
	/*                     if called.  Use  RCOND  in DGECO for a reliable */
	/*                     indication of singularity. */
	
	/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
	/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
	/* ***ROUTINES CALLED  DAXPY, DSCAL, IDAMAX */
	/* ***REVISION HISTORY  (YYMMDD) */
	/*   780814  DATE WRITTEN */
	/*   890831  Modified array declarations.  (WRB) */
	/*   890831  REVISION DATE from Version 3.2 */
	/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
	/*   900326  Removed duplicate information from DESCRIPTION section. */
	/*           (WRB) */
	/*   920501  Reformatted the REFERENCES section.  (WRB) */
	/* ***END PROLOGUE  DGEFA */
	
	
	/*     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING */
	
	/* ***FIRST EXECUTABLE STATEMENT  DGEFA */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    //a -= a_offset;
    //--ipvt;
	
    /* Function Body */
    *info = 0;
    nm1 = *PARAM_n - 1;
    if (nm1 < 1) {
		goto L70;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
		kp1 = k + 1;
		
		/*        FIND L = PIVOT INDEX */
		
		i__2 = *PARAM_n - k + 1;
		DGEFA_l = idamax_(&i__2, &a[k + k * a_dim1 -a_offset], 1, common) + k - 1;
		ipvt[k -1] = DGEFA_l;
		
		/*        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED */
		
		if (a[DGEFA_l + k * a_dim1 -a_offset] == 0.) {
			goto L40;
		}
		
		/*           INTERCHANGE IF NECESSARY */
		
		if (DGEFA_l == k) {
			goto L10;
		}
		t = a[DGEFA_l + k * a_dim1 -a_offset];
		a[DGEFA_l + k * a_dim1 -a_offset] = a[k + k * a_dim1 -a_offset];
		a[k + k * a_dim1 -a_offset] = t;
	L10:
		
		/*           COMPUTE MULTIPLIERS */
		
		t = -1. / a[k + k * a_dim1 -a_offset];
		i__2 = *PARAM_n - k;
		dscal_(&i__2, &t, &a[k + 1 + k * a_dim1 -a_offset], 1, common);
		
		/*           ROW ELIMINATION WITH COLUMN INDEXING */
		
		i__2 = *PARAM_n;
		for (j = kp1; j <= i__2; ++j) {
			t = a[DGEFA_l + j * a_dim1 -a_offset];
			if (DGEFA_l == k) {
				goto L20;
			}
			a[DGEFA_l + j * a_dim1 -a_offset] = a[k + j * a_dim1 -a_offset];
			a[k + j * a_dim1 -a_offset] = t;
		L20:
			i__3 = *PARAM_n - k;
			daxpy_(&i__3, &t, &a[k + 1 + k * a_dim1 -a_offset], 1, &a[k + 1 + j * a_dim1 -a_offset], 1, common);
			/* L30: */
		}
		goto L50;
	L40:
		*info = k;
	L50:
		/* L60: */
		;
    }
L70:
    ipvt[*PARAM_n -1] = *PARAM_n;
    if (a[*PARAM_n + *PARAM_n * a_dim1 -a_offset] == 0.) {
		*info = *PARAM_n;
    }
    return 0;
} /* dgefa_ */

/* DECK DGESL */
/* Subroutine */ 
__device__ int dgesl_(double  * __restrict__ a, int * __restrict__ lda, int * __restrict__ PARAM_n, int * __restrict__ ipvt, 
                      double  * __restrict__ b, int job, struct cuLsodaCommonBlock * __restrict__ common)
{
    /* System generated locals */
    int a_dim1 = 0;
	int a_offset = 0;
	int i__1 = 0;
	int i__2 = 0;
	
    /* Local variables */
    int k = 0;
	int DGESL_l = 0.;
    double  t = 0.;
    int kb = 0;
	int nm1 = 0;
	
	
	/* ***BEGIN PROLOGUE  DGESL */
	/* ***PURPOSE  Solve the real system A*X=B or TRANS(A)*X=B using the */
	/*            factors computed by DGECO or DGEFA. */
	/* ***CATEGORY  D2A1 */
	/* ***TYPE      DOUBLE PRECISION (SGESL-S, DGESL-D, CGESL-C) */
	/* ***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE */
	/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
	/* ***DESCRIPTION */
	
	/*     DGESL solves the double  precision system */
	/*     A * X = B  or  TRANS(A) * X = B */
	/*     using the factors computed by DGECO or DGEFA. */
	
	/*     On Entry */
	
	/*        A       DOUBLE PRECISION(LDA, N) */
	/*                the output from DGECO or DGEFA. */
	
	/*        LDA     int */
	/*                the leading dimension of the array  A . */
	
	/*        N       int */
	/*                the order of the matrix  A . */
	
	/*        IPVT    int(N) */
	/*                the pivot vector from DGECO or DGEFA. */
	
	/*        B       DOUBLE PRECISION(N) */
	/*                the right hand side vector. */
	
	/*        JOB     int */
	/*                = 0         to solve  A*X = B , */
	/*                = nonzero   to solve  TRANS(A)*X = B  where */
	/*                            TRANS(A)  is the transpose. */
	
	/*     On Return */
	
	/*        B       the solution vector  X . */
	
	/*     Error Condition */
	
	/*        A division by zero will occur if the input factor contains a */
	/*        zero on the diagonal.  Technically this indicates singularity */
	/*        but it is often caused by improper arguments or improper */
	/*        setting of LDA .  It will not occur if the subroutines are */
	/*        called correctly and if DGECO has set RCOND .GT. 0.0 */
	/*        or DGEFA has set INFO .EQ. 0 . */
	
	/*     To compute  INVERSE(A) * C  where  C  is a matrix */
	/*     with  P  columns */
	/*           CALL DGECO(A,LDA,N,IPVT,RCOND,Z) */
	/*           IF (RCOND is too small) GO TO ... */
	/*           DO 10 J = 1, P */
	/*              CALL DGESL(A,LDA,N,IPVT,C(1,J),0) */
	/*        10 CONTINUE */
	
	/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
	/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
	/* ***ROUTINES CALLED  DAXPY, DDOT */
	/* ***REVISION HISTORY  (YYMMDD) */
	/*   780814  DATE WRITTEN */
	/*   890831  Modified array declarations.  (WRB) */
	/*   890831  REVISION DATE from Version 3.2 */
	/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
	/*   900326  Removed duplicate information from DESCRIPTION section. */
	/*           (WRB) */
	/*   920501  Reformatted the REFERENCES section.  (WRB) */
	/* ***END PROLOGUE  DGESL */
	
	/* ***FIRST EXECUTABLE STATEMENT  DGESL */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    //a -= a_offset;
    //--ipvt;
    //--b;
	
    /* Function Body */
    nm1 = *PARAM_n - 1;
    if (job != 0) {
		goto L50;
    }
	
	/*        JOB = 0 , SOLVE  A * X = B */
	/*        FIRST SOLVE  L*Y = B */
	
    if (nm1 < 1) {
		goto L30;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
		DGESL_l = ipvt[k -1];
		t = b[DGESL_l -1];
		if (DGESL_l == k) {
			goto L10;
		}
		b[DGESL_l -1] = b[k -1];
		b[k -1] = t;
	L10:
		i__2 = *PARAM_n - k;
		daxpy_(&i__2, &t, &a[k + 1 + k * a_dim1 -a_offset], 1, &b[k + 1 -1], 1, common);
		/* L20: */
    }
L30:
	
	/*        NOW SOLVE  U*X = Y */
	
    i__1 = *PARAM_n;
    for (kb = 1; kb <= i__1; ++kb) {
		k = *PARAM_n + 1 - kb;
		b[k -1] /= a[k + k * a_dim1 -a_offset];
		t = -b[k -1];
		i__2 = k - 1;
		daxpy_(&i__2, &t, &a[k * a_dim1 + 1 -a_offset], 1, b, 1, common);
		/* L40: */
    }
    goto L100;
L50:
	
	/*        JOB = NONZERO, SOLVE  TRANS(A) * X = B */
	/*        FIRST SOLVE  TRANS(U)*Y = B */
	
    i__1 = *PARAM_n;
    for (k = 1; k <= i__1; ++k) {
		i__2 = k - 1;
		t = ddot_(&i__2, &a[k * a_dim1 + 1 -a_offset], 1, b, 1, common);
		b[k -1] = (b[k -1] - t) / a[k + k * a_dim1 -a_offset];
		/* L60: */
    }
	
	/*        NOW SOLVE TRANS(L)*X = Y */
	
    if (nm1 < 1) {
		goto L90;
    }
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
		k = *PARAM_n - kb;
		i__2 = *PARAM_n - k;
		b[k -1] += ddot_(&i__2, &a[k + 1 + k * a_dim1 -a_offset], 1, &b[k + 1 -1], 1, common);
		DGESL_l = ipvt[k -1];
		if (DGESL_l == k) {
			goto L70;
		}
		t = b[DGESL_l -1];
		b[DGESL_l -1] = b[k -1];
		b[k -1] = t;
	L70:
		/* L80: */
		;
    }
L90:
L100:
    return 0;
} /* dgesl_ */

/* DECK DGBFA */
/* Subroutine */ 
__device__ int dgbfa_(double  * __restrict__ abd, int * __restrict__ lda, int * __restrict__ PARAM_n, 
                      int * __restrict__ ml, int * __restrict__ mu, 
                      int * __restrict__ ipvt, int * __restrict__ info, struct cuLsodaCommonBlock * __restrict__ common)
{
    /* System generated locals */
    int abd_dim1 = 0;
	int abd_offset = 0;
	int i__1 = 0;
	int i__2 = 0;
	int i__3 = 0;
	int i__4 = 0;
	
    /* Local variables */
    int i__ = 0;
	int j = 0;
	int k = 0;
	int DGBFA_l = 0;
	int m = 0;
    double  t = 0.;
    int i0 = 0;
	int j0 = 0;
	int j1 = 0;
	int lm = 0;
	int mm = 0;
	int ju = 0;
	int jz = 0;
	int kp1 = 0;
	int nm1 = 0;
	
	/* ***BEGIN PROLOGUE  DGBFA */
	/* ***PURPOSE  Factor a band matrix using Gaussian elimination. */
	/* ***CATEGORY  D2A2 */
	/* ***TYPE      DOUBLE PRECISION (SGBFA-S, DGBFA-D, CGBFA-C) */
	/* ***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION */
	/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
	/* ***DESCRIPTION */
	
	/*     DGBFA factors a double  precision band matrix by elimination. */
	
	/*     DGBFA is usually called by DGBCO, but it can be called */
	/*     directly with a saving in time if  RCOND  is not needed. */
	
	/*     On Entry */
	
	/*        ABD     DOUBLE PRECISION(LDA, N) */
	/*                contains the matrix in band storage.  The columns */
	/*                of the matrix are stored in the columns of  ABD  and */
	/*                the diagonals of the matrix are stored in rows */
	/*                ML+1 through 2*ML+MU+1 of  ABD . */
	/*                See the comments below for details. */
	
	/*        LDA     int */
	/*                the leading dimension of the array  ABD . */
	/*                LDA must be .GE. 2*ML + MU + 1 . */
	
	/*        N       int */
	/*                the order of the original matrix. */
	
	/*        ML      int */
	/*                number of diagonals below the main diagonal. */
	/*                0 .LE. ML .LT.  N . */
	
	/*        MU      int */
	/*                number of diagonals above the main diagonal. */
	/*                0 .LE. MU .LT.  N . */
	/*                More efficient if  ML .LE. MU . */
	/*     On Return */
	
	/*        ABD     an upper triangular matrix in band storage and */
	/*                the multipliers which were used to obtain it. */
	/*                The factorization can be written  A = L*U  where */
	/*                L  is a product of permutation and unit lower */
	/*                triangular matrices and  U  is upper triangular. */
	
	/*        IPVT    int(N) */
	/*                an int vector of pivot indices. */
	
	/*        INFO    int */
	/*                = 0  normal value. */
	/*                = K  if  U(K,K) .EQ. 0.0 .  This is not an error */
	/*                     condition for this subroutine, but it does */
	/*                     indicate that DGBSL will divide by zero if */
	/*                     called.  Use  RCOND  in DGBCO for a reliable */
	/*                     indication of singularity. */
	
	/*     Band Storage */
	
	/*           If  A  is a band matrix, the following program segment */
	/*           will set up the input. */
	
	/*                   ML = (band width below the diagonal) */
	/*                   MU = (band width above the diagonal) */
	/*                   M = ML + MU + 1 */
	/*                   DO 20 J = 1, N */
	/*                      I1 = MAX(1, J-MU) */
	/*                      I2 = MIN(N, J+ML) */
	/*                      DO 10 I = I1, I2 */
	/*                         K = I - J + M */
	/*                         ABD(K,J) = A(I,J) */
	/*                10    CONTINUE */
	/*                20 CONTINUE */
	
	/*           This uses rows  ML+1  through  2*ML+MU+1  of  ABD . */
	/*           In addition, the first  ML  rows in  ABD  are used for */
	/*           elements generated during the triangularization. */
	/*           The total number of rows needed in  ABD  is  2*ML+MU+1 . */
	/*           The  ML+MU by ML+MU  upper left triangle and the */
	/*           ML by ML  lower right triangle are not referenced. */
	
	/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
	/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
	/* ***ROUTINES CALLED  DAXPY, DSCAL, IDAMAX */
	/* ***REVISION HISTORY  (YYMMDD) */
	/*   780814  DATE WRITTEN */
	/*   890531  Changed all specific intrinsics to generic.  (WRB) */
	/*   890831  Modified array declarations.  (WRB) */
	/*   890831  REVISION DATE from Version 3.2 */
	/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
	/*   900326  Removed duplicate information from DESCRIPTION section. */
	/*           (WRB) */
	/*   920501  Reformatted the REFERENCES section.  (WRB) */
	/* ***END PROLOGUE  DGBFA */
	
	
	/* ***FIRST EXECUTABLE STATEMENT  DGBFA */
    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    //abd -= abd_offset;
    //--ipvt;
	
    /* Function Body */
    m = *ml + *mu + 1;
    *info = 0;
	
	/*     ZERO INITIAL FILL-IN COLUMNS */
	
    j0 = *mu + 2;
    j1 = min(*PARAM_n,m) - 1;
    if (j1 < j0) {
		goto L30;
    }
    i__1 = j1;
    for (jz = j0; jz <= i__1; ++jz) {
		i0 = m + 1 - jz;
		i__2 = *ml;
		for (i__ = i0; i__ <= i__2; ++i__) {
			abd[i__ + jz * abd_dim1 -abd_offset] = 0.;
			/* L10: */
		}
		/* L20: */
    }
L30:
    jz = j1;
    ju = 0;
	
	/*     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING */
	
    nm1 = *PARAM_n - 1;
    if (nm1 < 1) {
		goto L130;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
		kp1 = k + 1;
		
		/*        ZERO NEXT FILL-IN COLUMN */
		
		++jz;
		if (jz > *PARAM_n) {
			goto L50;
		}
		if (*ml < 1) {
			goto L50;
		}
		i__2 = *ml;
		for (i__ = 1; i__ <= i__2; ++i__) {
			abd[i__ + jz * abd_dim1 -abd_offset] = 0.;
			/* L40: */
		}
	L50:
		
		/*        FIND L = PIVOT INDEX */
		
		/* Computing MIN */
		i__2 = *ml, i__3 = *PARAM_n - k;
		lm = min(i__2,i__3);
		i__2 = lm + 1;
		DGBFA_l = idamax_(&i__2, &abd[m + k * abd_dim1 -abd_offset], 1, common) + m - 1;
		ipvt[k -1] = DGBFA_l + k - m;
		
		/*        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED */
		
		if (abd[DGBFA_l + k * abd_dim1 -abd_offset] == 0.) {
			goto L100;
		}
		
		/*           INTERCHANGE IF NECESSARY */
		
		if (DGBFA_l == m) {
			goto L60;
		}
		t = abd[DGBFA_l + k * abd_dim1 -abd_offset];
		abd[DGBFA_l + k * abd_dim1 -abd_offset] = abd[m + k * abd_dim1 -abd_offset];
		abd[m + k * abd_dim1 -abd_offset] = t;
	L60:
		
		/*           COMPUTE MULTIPLIERS */
		
		t = -1. / abd[m + k * abd_dim1 -abd_offset];
		dscal_(&lm, &t, &abd[m + 1 + k * abd_dim1 -abd_offset], 1, common);
		
		/*           ROW ELIMINATION WITH COLUMN INDEXING */
		
		/* Computing MIN */
		/* Computing MAX */
		i__3 = ju, i__4 = *mu + ipvt[k -1];
		i__2 = max(i__3,i__4);
		ju = min(i__2,*PARAM_n);
		mm = m;
		if (ju < kp1) {
			goto L90;
		}
		i__2 = ju;
		for (j = kp1; j <= i__2; ++j) {
			--DGBFA_l;
			--mm;
			t = abd[DGBFA_l + j * abd_dim1 -abd_offset];
			if (DGBFA_l == mm) {
				goto L70;
			}
			abd[DGBFA_l + j * abd_dim1 -abd_offset] = abd[mm + j * abd_dim1 -abd_offset];
			abd[mm + j * abd_dim1 -abd_offset] = t;
		L70:
			daxpy_(&lm, &t, &abd[m + 1 + k * abd_dim1 -abd_offset], 1, &abd[mm + 1 + j * abd_dim1 -abd_offset], 1, common);
			/* L80: */
		}
	L90:
		goto L110;
	L100:
		*info = k;
	L110:
		/* L120: */
		;
    }
L130:
    ipvt[*PARAM_n -1] = *PARAM_n;
    if (abd[m + *PARAM_n * abd_dim1 -abd_offset] == 0.) {
		*info = *PARAM_n;
    }
    return 0;
} /* dgbfa_ */

/* DECK DGBSL */
/* Subroutine */ 
__device__ int dgbsl_(double  * __restrict__ abd, int * __restrict__ lda, int * __restrict__ PARAM_n, 
                      int * __restrict__ ml, int * __restrict__ mu, int * __restrict__ ipvt, double  * __restrict__ b, int job, struct cuLsodaCommonBlock *common)
{
    /* System generated locals */
    int abd_dim1 = 0;
	int abd_offset = 0;
	int i__1 = 0;
	int i__2 = 0;
	int i__3 = 0;

	
    /* Local variables */
    int k = 0;
	int DGBSL_l = 0;
	int m = 0;
    double  t = 0.;
    int kb = 0;
	int la = 0;
	int lb = 0;
	int lm = 0;
	int nm1 = 0;
	
	
	/* ***BEGIN PROLOGUE  DGBSL */
	/* ***PURPOSE  Solve the real band system A*X=B or TRANS(A)*X=B using */
	/*            the factors computed by DGBCO or DGBFA. */
	/* ***CATEGORY  D2A2 */
	/* ***TYPE      DOUBLE PRECISION (SGBSL-S, DGBSL-D, CGBSL-C) */
	/* ***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE */
	/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
	/* ***DESCRIPTION */
	
	/*     DGBSL solves the double  precision band system */
	/*     A * X = B  or  TRANS(A) * X = B */
	/*     using the factors computed by DGBCO or DGBFA. */
	
	/*     On Entry */
	
	/*        ABD     DOUBLE PRECISION(LDA, N) */
	/*                the output from DGBCO or DGBFA. */
	
	/*        LDA     int */
	/*                the leading dimension of the array  ABD . */
	
	/*        N       int */
	/*                the order of the original matrix. */
	
	/*        ML      int */
	/*                number of diagonals below the main diagonal. */
	
	/*        MU      int */
	/*                number of diagonals above the main diagonal. */
	
	/*        IPVT    int(N) */
	/*                the pivot vector from DGBCO or DGBFA. */
	
	/*        B       DOUBLE PRECISION(N) */
	/*                the right hand side vector. */
	
	/*        JOB     int */
	/*                = 0         to solve  A*X = B , */
	/*                = nonzero   to solve  TRANS(A)*X = B , where */
	/*                            TRANS(A)  is the transpose. */
	
	/*     On Return */
	
	/*        B       the solution vector  X . */
	
	/*     Error Condition */
	
	/*        A division by zero will occur if the input factor contains a */
	/*        zero on the diagonal.  Technically this indicates singularity */
	/*        but it is often caused by improper arguments or improper */
	/*        setting of LDA .  It will not occur if the subroutines are */
	/*        called correctly and if DGBCO has set RCOND .GT. 0.0 */
	/*        or DGBFA has set INFO .EQ. 0 . */
	
	/*     To compute  INVERSE(A) * C  where  C  is a matrix */
	/*     with  P  columns */
	/*           CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z) */
	/*           IF (RCOND is too small) GO TO ... */
	/*           DO 10 J = 1, P */
	/*              CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0) */
	/*        10 CONTINUE */
	
	/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
	/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
	/* ***ROUTINES CALLED  DAXPY, DDOT */
	/* ***REVISION HISTORY  (YYMMDD) */
	/*   780814  DATE WRITTEN */
	/*   890531  Changed all specific intrinsics to generic.  (WRB) */
	/*   890831  Modified array declarations.  (WRB) */
	/*   890831  REVISION DATE from Version 3.2 */
	/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
	/*   900326  Removed duplicate information from DESCRIPTION section. */
	/*           (WRB) */
	/*   920501  Reformatted the REFERENCES section.  (WRB) */
	/* ***END PROLOGUE  DGBSL */
	
	/* ***FIRST EXECUTABLE STATEMENT  DGBSL */
    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    //abd -= abd_offset;
    //--ipvt;
    //--b;
	
    /* Function Body */
    m = *mu + *ml + 1;
    nm1 = *PARAM_n - 1;
    if (job != 0) {
		goto L50;
    }
	
	/*        JOB = 0 , SOLVE  A * X = B */
	/*        FIRST SOLVE L*Y = B */
	
    if (*ml == 0) {
		goto L30;
    }
    if (nm1 < 1) {
		goto L30;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
		/* Computing MIN */
		i__2 = *ml, i__3 = *PARAM_n - k;
		lm = min(i__2,i__3);
		DGBSL_l = ipvt[k -1];
		t = b[DGBSL_l -1];
		if (DGBSL_l == k) {
			goto L10;
		}
		b[DGBSL_l -1] = b[k -1];
		b[k -1] = t;
	L10:
		daxpy_(&lm, &t, &abd[m + 1 + k * abd_dim1 -abd_offset], 1, &b[k + 1 -1], 1, common);
		/* L20: */
    }
L30:
	
	/*        NOW SOLVE  U*X = Y */
	
    i__1 = *PARAM_n;
    for (kb = 1; kb <= i__1; ++kb) {
		k = *PARAM_n + 1 - kb;
		b[k -1] /= abd[m + k * abd_dim1 -abd_offset];
		lm = min(k,m) - 1;
		la = m - lm;
		lb = k - lm;
		t = -b[k -1];
		daxpy_(&lm, &t, &abd[la + k * abd_dim1 -abd_offset], 1, &b[lb -1], 1, common);
		/* L40: */
    }
    goto L100;
L50:
	
	/*        JOB = NONZERO, SOLVE  TRANS(A) * X = B */
	/*        FIRST SOLVE  TRANS(U)*Y = B */
	
    i__1 = *PARAM_n;
    for (k = 1; k <= i__1; ++k) {
		lm = min(k,m) - 1;
		la = m - lm;
		lb = k - lm;
		t = ddot_(&lm, &abd[la + k * abd_dim1 -abd_offset], 1, &b[lb -1], 1, common);
		b[k -1] = (b[k -1] - t) / abd[m + k * abd_dim1 -abd_offset];
		/* L60: */
    }
	
	/*        NOW SOLVE TRANS(L)*X = Y */
	
    if (*ml == 0) {
		goto L90;
    }
    if (nm1 < 1) {
		goto L90;
    }
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
		k = *PARAM_n - kb;
		/* Computing MIN */
		i__2 = *ml, i__3 = *PARAM_n - k;
		lm = min(i__2,i__3);
		b[k -1] += ddot_(&lm, &abd[m + 1 + k * abd_dim1 -abd_offset], 1, &b[k + 1 -1], 1, common);
		DGBSL_l = ipvt[k -1];
		if (DGBSL_l == k) {
			goto L70;
		}
		t = b[DGBSL_l -1];
		b[DGBSL_l -1] = b[k -1];
		b[k -1] = t;
	L70:
		/* L80: */
		;
    }
L90:
L100:
    return 0;
} /* dgbsl_ */

/* DECK DUMACH */
__device__ double  dumach_(struct cuLsodaCommonBlock * __restrict__ common)
{
    /* System generated locals */
    double  ret_val = 0.;
	
    /* Local variables */
	double  u = 0.;
	double  comp = 0.;
	
	/* ***BEGIN PROLOGUE  DUMACH */
	/* ***PURPOSE  Compute the unit roundoff of the machine. */
	/* ***CATEGORY  R1 */
	/* ***TYPE      DOUBLE PRECISION (RUMACH-S, DUMACH-D) */
	/* ***KEYWORDS  MACHINE CONSTANTS */
	/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
	/* ***DESCRIPTION */
	/* *Usage: */
	/*        DOUBLE PRECISION  A, DUMACH */
	/*        A = DUMACH() */
	
	/* *Function Return Values: */
	/*     A : the unit roundoff of the machine. */
	
	/* *Description: */
	/*     The unit roundoff is defined as the smallest positive machine */
	/*     number u such that  1.0 + u .ne. 1.0.  This is computed by DUMACH */
	/*     in a machine-independent manner. */
	
	/* ***REFERENCES  (NONE) */
	/* ***ROUTINES CALLED  DUMSUM */
	/* ***REVISION HISTORY  (YYYYMMDD) */
	/*   19930216  DATE WRITTEN */
	/*   19930818  Added SLATEC-format prologue.  (FNF) */
	/*   20030707  Added DUMSUM to force normal storage of COMP.  (ACH) */
	/* ***END PROLOGUE  DUMACH */
	
	/* ***FIRST EXECUTABLE STATEMENT  DUMACH */
    u = 1.;
L10:
    u *= .5;
    dumsum_(1., u, &comp, common);
    if (comp != 1.) {
		goto L10;
    }
    ret_val = u * 2.;
    return ret_val;
	/* ----------------------- End of Function DUMACH ------------------------ */
} /* dumach_ */

/* DECK XSETF */
/* Subroutine */ 
//__device__ int xsetf_(int *mflag, struct cuLsodaCommonBlock *common)
//{
//     int junk;
	
	
	/* ***BEGIN PROLOGUE  XSETF */
	/* ***PURPOSE  Reset the error print control flag. */
	/* ***CATEGORY  R3A */
	/* ***TYPE      ALL (XSETF-A) */
	/* ***KEYWORDS  ERROR CONTROL */
	/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
	/* ***DESCRIPTION */
	
	/*   XSETF sets the error print control flag to MFLAG: */
	/*      MFLAG=1 means print all messages (the default). */
	/*      MFLAG=0 means no printing. */
	
	/* ***SEE ALSO  XERRWD, XERRWV */
	/* ***REFERENCES  (NONE) */
	/* ***ROUTINES CALLED  IXSAV */
	/* ***REVISION HISTORY  (YYMMDD) */
	/*   921118  DATE WRITTEN */
	/*   930329  Added SLATEC format prologue. (FNF) */
	/*   930407  Corrected SEE ALSO section. (FNF) */
	/*   930922  Made user-callable, and other cosmetic changes. (FNF) */
	/* ***END PROLOGUE  XSETF */
	
	/* Subroutines called by XSETF.. None */
	/* Function routine called by XSETF.. IXSAV */
	/* ----------------------------------------------------------------------- */
	/* **End */
	
	/* ***FIRST EXECUTABLE STATEMENT  XSETF */
//    if (*mflag == 0 || *mflag == 1) {
//		junk = ixsav_(2, mflag, 1, common);
//    }
//    return 0;
	/* ----------------------- End of Subroutine XSETF ----------------------- */
//} /* xsetf_ */

/* DECK XSETUN */
/* Subroutine */ 
//__device__ int xsetun_(int *lun, struct cuLsodaCommonBlock *common)
//{
//     int junk;
	
	
	/* ***BEGIN PROLOGUE  XSETUN */
	/* ***PURPOSE  Reset the logical unit number for error messages. */
	/* ***CATEGORY  R3B */
	/* ***TYPE      ALL (XSETUN-A) */
	/* ***KEYWORDS  ERROR CONTROL */
	/* ***DESCRIPTION */
	
	/*   XSETUN sets the logical unit number for error messages to LUN. */
	
	/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
	/* ***SEE ALSO  XERRWD, XERRWV */
	/* ***REFERENCES  (NONE) */
	/* ***ROUTINES CALLED  IXSAV */
	/* ***REVISION HISTORY  (YYMMDD) */
	/*   921118  DATE WRITTEN */
	/*   930329  Added SLATEC format prologue. (FNF) */
	/*   930407  Corrected SEE ALSO section. (FNF) */
	/*   930922  Made user-callable, and other cosmetic changes. (FNF) */
	/* ***END PROLOGUE  XSETUN */
	
	/* Subroutines called by XSETUN.. None */
	/* Function routine called by XSETUN.. IXSAV */
	/* ----------------------------------------------------------------------- */
	/* **End */
	
	/* ***FIRST EXECUTABLE STATEMENT  XSETUN */
//    if (*lun > 0) {
//		junk = ixsav_(1, lun, 1, common);
//    }
//    return 0;
	/* ----------------------- End of Subroutine XSETUN ---------------------- */
//} /* xsetun_ */

/* DECK IXSAV */
__device__ int ixsav_(int ipar, int * __restrict__ ivalue, int iset, struct cuLsodaCommonBlock * __restrict__ common)
{
    /* Initialized data */
	
     int lunit = -1;
     int mesflg = 1;
	
    /* System generated locals */
    int ret_val = 0;
	
    /* Local variables */
	
	
	/* ***BEGIN PROLOGUE  IXSAV */
	/* ***SUBSIDIARY */
	/* ***PURPOSE  Save and recall error message control parameters. */
	/* ***CATEGORY  R3C */
	/* ***TYPE      ALL (IXSAV-A) */
	/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
	/* ***DESCRIPTION */
	
	/*  IXSAV saves and recalls one of two error message parameters: */
	/*    LUNIT, the logical unit number to which messages are printed, and */
	/*    MESFLG, the message print flag. */
	/*  This is a modification of the SLATEC library routine J4SAVE. */
	
	/*  Saved local variables.. */
	/*   LUNIT  = Logical unit number for messages.  The default is obtained */
	/*            by a call to IUMACH (may be machine-dependent). */
	/*   MESFLG = Print control flag.. */
	/*            1 means print all messages (the default). */
	/*            0 means no printing. */
	
	/*  On input.. */
	/*    IPAR   = Parameter indicator (1 for LUNIT, 2 for MESFLG). */
	/*    IVALUE = The value to be set for the parameter, if ISET = .TRUE. */
	/*    ISET   = Logical flag to indicate whether to read or write. */
	/*             If ISET = .TRUE., the parameter will be given */
	/*             the value IVALUE.  If ISET = .FALSE., the parameter */
	/*             will be unchanged, and IVALUE is a dummy argument. */
	
	/*  On return.. */
	/*    IXSAV = The (old) value of the parameter. */
	
	/* ***SEE ALSO  XERRWD, XERRWV */
	/* ***ROUTINES CALLED  IUMACH */
	/* ***REVISION HISTORY  (YYMMDD) */
	/*   921118  DATE WRITTEN */
	/*   930329  Modified prologue to SLATEC format. (FNF) */
	/*   930915  Added IUMACH call to get default output unit.  (ACH) */
	/*   930922  Minor cosmetic changes. (FNF) */
	/*   010425  Type declaration for IUMACH added. (ACH) */
	/* ***END PROLOGUE  IXSAV */
	
	/* Subroutines called by IXSAV.. None */
	/* Function routine called by IXSAV.. IUMACH */
	/* ----------------------------------------------------------------------- */
	/* **End */
	/* ----------------------------------------------------------------------- */
	/* ----------------------------------------------------------------------- */
	/* The following Fortran-77 declaration is to cause the values of the */
	/* listed (local) variables to be saved between calls to this routine. */
	/* ----------------------------------------------------------------------- */
	
	/* ***FIRST EXECUTABLE STATEMENT  IXSAV */
    if (ipar == 1) {
		if (lunit == -1) {
			lunit = 6;
		}
		ret_val = lunit;
		if (iset) {
			lunit = *ivalue;
		}
    }
	
    if (ipar == 2) {
		ret_val = mesflg;
		if (iset) { 
			mesflg = *ivalue;
		}
    }
	
    return ret_val;
	/* ----------------------- End of Function IXSAV ------------------------- */
} /* ixsav_ */


/* DECK IDAMAX */
__device__ int idamax_(int *PARAM_n, double  *dx, int incx, struct cuLsodaCommonBlock *common)
{
    /* System generated locals */
    int ret_val = 0;
	int i__1 = 0;
    double  d__1 = 0;
	
    /* Local variables */
    int i__ = 0;
	int ix = 0;
    double  dmax__ = 0.;
	double  xmag = 0.;
	
	/* ***BEGIN PROLOGUE  IDAMAX */
	/* ***PURPOSE  Find the smallest index of that component of a vector */
	/*            having the maximum magnitude. */
	/* ***CATEGORY  D1A2 */
	/* ***TYPE      DOUBLE PRECISION (ISAMAX-S, IDAMAX-D, ICAMAX-C) */
	/* ***KEYWORDS  BLAS, LINEAR ALGEBRA, MAXIMUM COMPONENT, VECTOR */
	/* ***AUTHOR  Lawson, C. L., (JPL) */
	/*           Hanson, R. J., (SNLA) */
	/*           Kincaid, D. R., (U. of Texas) */
	/*           Krogh, F. T., (JPL) */
	/* ***DESCRIPTION */
	
	/*                B L A S  Subprogram */
	/*    Description of Parameters */
	
	/*     --Input-- */
	/*        N  number of elements in input vector(s) */
	/*       DX  double  precision vector with N elements */
	/*     INCX  storage spacing between elements of DX */
	
	/*     --Output-- */
	/*   IDAMAX  smallest index (zero if N .LE. 0) */
	
	/*     Find smallest index of maximum magnitude of double  precision DX. */
	/*     IDAMAX = first I, I = 1 to N, to maximize ABS(DX(IX+(I-1)*INCX)), */
	/*     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX. */
	
	/* ***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. */
	/*                 Krogh, Basic linear algebra subprograms for Fortran */
	/*                 usage, Algorithm No. 539, Transactions on Mathematical */
	/*                 Software 5, 3 (September 1979), pp. 308-323. */
	/* ***ROUTINES CALLED  (NONE) */
	/* ***REVISION HISTORY  (YYMMDD) */
	/*   791001  DATE WRITTEN */
	/*   890531  Changed all specific intrinsics to generic.  (WRB) */
	/*   890531  REVISION DATE from Version 3.2 */
	/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
	/*   900821  Modified to correct problem with a negative increment. */
	/*           (WRB) */
	/*   920501  Reformatted the REFERENCES section.  (WRB) */
	/* ***END PROLOGUE  IDAMAX */
	/* ***FIRST EXECUTABLE STATEMENT  IDAMAX */
   
	 /* Parameter adjustments */
    //--dx;
	
    /* Function Body */
    ret_val = 0;
    if (*PARAM_n <= 0) {
		return ret_val;
    }
    ret_val = 1;
    if (*PARAM_n == 1) {
		return ret_val;
    }
	
    if (incx == 1) {
		goto L20;
    }
	
	/*     Code for increments not equal to 1. */
	
    ix = 1;
    if (incx < 0) {
		ix = (-(*PARAM_n) + 1) * incx + 1;
    }
    dmax__ = (d__1 = dx[ix -1], fabs(d__1));
    ix += incx;
    i__1 = *PARAM_n;
    for (i__ = 2; i__ <= i__1; ++i__) {
		xmag = (d__1 = dx[ix -1], fabs(d__1));
		if (xmag > dmax__) {
			ret_val = i__;
			dmax__ = xmag;
		}
		ix += incx;
		/* L10: */
    }
    return ret_val;
	
	/*     Code for increments equal to 1. */
	
L20:
    dmax__ = fabs(dx[0]);
    i__1 = *PARAM_n;
    for (i__ = 2; i__ <= i__1; ++i__) {
		xmag = (d__1 = dx[i__ -1], fabs(d__1));
		if (xmag > dmax__) {
			ret_val = i__;
			dmax__ = xmag;
		}
		/* L30: */
    }
    return ret_val;
} /* idamax_ */

/* DECK DAXPY */
/* Subroutine */
__device__ int daxpy_(int * __restrict__ PARAM_n, double  * __restrict__ da, double  * __restrict__ dx, 
                      int incx, double  * __restrict__ dy, int incy, struct cuLsodaCommonBlock * __restrict__ common)
{
    /* System generated locals */
    int i__1 = 0;
	int i__2 = 0;
	
    /* Local variables */
	int i__ = 0;
	int m = 0;
	int ix = 0;
	int iy = 0;
	int ns = 0;
	int mp1 = 0;
	
	/* ***BEGIN PROLOGUE  DAXPY */
	/* ***PURPOSE  Compute a constant times a vector plus a vector. */
	/* ***CATEGORY  D1A7 */
	/* ***TYPE      DOUBLE PRECISION (SAXPY-S, DAXPY-D, CAXPY-C) */
	/* ***KEYWORDS  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR */
	/* ***AUTHOR  Lawson, C. L., (JPL) */
	/*           Hanson, R. J., (SNLA) */
	/*           Kincaid, D. R., (U. of Texas) */
	/*           Krogh, F. T., (JPL) */
	/* ***DESCRIPTION */
	
	/*                B L A S  Subprogram */
	/*    Description of Parameters */
	
	/*     --Input-- */
	/*        N  number of elements in input vector(s) */
	/*       DA  double  precision scalar multiplier */
	/*       DX  double  precision vector with N elements */
	/*     INCX  storage spacing between elements of DX */
	/*       DY  double  precision vector with N elements */
	/*     INCY  storage spacing between elements of DY */
	
	/*     --Output-- */
	/*       DY  double  precision result (unchanged if N .LE. 0) */
	
	/*     Overwrite double  precision DY with double  precision DA*DX + DY. */
	/*     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) + */
	/*       DY(LY+I*INCY), */
	/*     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is */
	/*     defined in a similar way using INCY. */
	
	/* ***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. */
	/*                 Krogh, Basic linear algebra subprograms for Fortran */
	/*                 usage, Algorithm No. 539, Transactions on Mathematical */
	/*                 Software 5, 3 (September 1979), pp. 308-323. */
	/* ***ROUTINES CALLED  (NONE) */
	/* ***REVISION HISTORY  (YYMMDD) */
	/*   791001  DATE WRITTEN */
	/*   890831  Modified array declarations.  (WRB) */
	/*   890831  REVISION DATE from Version 3.2 */
	/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
	/*   920310  Corrected definition of LX in DESCRIPTION.  (WRB) */
	/*   920501  Reformatted the REFERENCES section.  (WRB) */
	/* ***END PROLOGUE  DAXPY */
	/* ***FIRST EXECUTABLE STATEMENT  DAXPY */
    /* Parameter adjustments */
    //--dy;
    //--dx;
	
    /* Function Body */
    if (*PARAM_n <= 0 || *da == 0.) {
		return 0;
    }
    if (incx == incy) {
		if ((i__1 = incx - 1) < 0) {
			goto L5;
		} else if (i__1 == 0) {
			goto L20;
		} else {
			goto L60;
		}
    }
	
	/*     Code for unequal or nonpositive increments. */
	
L5:
    ix = 1;
    iy = 1;
    if (incx < 0) {
		ix = (-(*PARAM_n) + 1) * incx + 1;
    }
    if (incy < 0) {
		iy = (-(*PARAM_n) + 1) * incy + 1;
    }
    i__1 = *PARAM_n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		dy[iy -1] += *da * dx[ix -1];
		ix += incx;
		iy += incy;
		/* L10: */
    }
    return 0;
	
	/*     Code for both increments equal to 1. */
	
	/*     Clean-up loop so remaining vector length is a multiple of 4. */
	
L20:
    m = *PARAM_n % 4;
    if (m == 0) {
		goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
		dy[i__ -1] += *da * dx[i__ -1];
		/* L30: */
    }
    if (*PARAM_n < 4) {
		return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *PARAM_n;
    for (i__ = mp1; i__ <= i__1; i__ += 4) {
		dy[i__ -1] += *da * dx[i__ -1];
		dy[i__ + 1 -1] += *da * dx[i__ + 1 -1];
		dy[i__ + 2 -1] += *da * dx[i__ + 2 -1];
		dy[i__ + 3 -1] += *da * dx[i__ + 3 -1];
		/* L50: */
    }
    return 0;
	
	/*     Code for equal, positive, non-unit increments. */
	
L60:
    ns = *PARAM_n * incx;
    i__1 = ns;
    i__2 = incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
		dy[i__ -1] = *da * dx[i__ -1] + dy[i__ -1];
		/* L70: */
    }
    return 0;
} /* daxpy_ */

/* Subroutine */ 
__device__ int dumsum_(double  a, double  b, double  * __restrict__ c__, struct cuLsodaCommonBlock * __restrict__ common)
{
	/*     Routine to force normal storing of A + B, for DUMACH. */
    *c__ = a + b;
    return 0;
} /* dumsum_ */

/* DECK DSCAL */
/* Subroutine */ 
__device__ int dscal_(int * __restrict__ PARAM_n, double  * __restrict__ da, double  * __restrict__dx, int incx, struct cuLsodaCommonBlock * __restrict__ common)
{
    /* System generated locals */
    int i__1 = 0;
	
    /* Local variables */
	int i__ = 0;
	int m = 0;
	int ix = 0;
	int mp1 = 0;
	
	/* ***BEGIN PROLOGUE  DSCAL */
	/* ***PURPOSE  Multiply a vector by a constant. */
	/* ***CATEGORY  D1A6 */
	/* ***TYPE      DOUBLE PRECISION (SSCAL-S, DSCAL-D, CSCAL-C) */
	/* ***KEYWORDS  BLAS, LINEAR ALGEBRA, SCALE, VECTOR */
	/* ***AUTHOR  Lawson, C. L., (JPL) */
	/*           Hanson, R. J., (SNLA) */
	/*           Kincaid, D. R., (U. of Texas) */
	/*           Krogh, F. T., (JPL) */
	/* ***DESCRIPTION */
	
	/*                B L A S  Subprogram */
	/*    Description of Parameters */
	
	/*     --Input-- */
	/*        N  number of elements in input vector(s) */
	/*       DA  double  precision scale factor */
	/*       DX  double  precision vector with N elements */
	/*     INCX  storage spacing between elements of DX */
	
	/*     --Output-- */
	/*       DX  double  precision result (unchanged if N.LE.0) */
	
	/*     Replace double  precision DX by double  precision DA*DX. */
	/*     For I = 0 to N-1, replace DX(IX+I*INCX) with  DA * DX(IX+I*INCX), */
	/*     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX. */
	
	/* ***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. */
	/*                 Krogh, Basic linear algebra subprograms for Fortran */
	/*                 usage, Algorithm No. 539, Transactions on Mathematical */
	/*                 Software 5, 3 (September 1979), pp. 308-323. */
	/* ***ROUTINES CALLED  (NONE) */
	/* ***REVISION HISTORY  (YYMMDD) */
	/*   791001  DATE WRITTEN */
	/*   890831  Modified array declarations.  (WRB) */
	/*   890831  REVISION DATE from Version 3.2 */
	/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
	/*   900821  Modified to correct problem with a negative increment. */
	/*           (WRB) */
	/*   920501  Reformatted the REFERENCES section.  (WRB) */
	/* ***END PROLOGUE  DSCAL */
	/* ***FIRST EXECUTABLE STATEMENT  DSCAL */
   
	 /* Parameter adjustments */
    //--dx;
	
    /* Function Body */
    if (*PARAM_n <= 0) {
		return 0;
    }
    if (incx == 1) {
		goto L20;
    }
	
	/*     Code for increment not equal to 1. */
	
    ix = 1;
    if (incx < 0) {
		ix = (-(*PARAM_n) + 1) * incx + 1;
    }
    i__1 = *PARAM_n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		dx[ix -1] = *da * dx[ix -1];
		ix += incx;
		/* L10: */
    }
    return 0;
	
	/*     Code for increment equal to 1. */
	
	/*     Clean-up loop so remaining vector length is a multiple of 5. */
	
L20:
    m = *PARAM_n % 5;
    if (m == 0) {
		goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
		dx[i__ -1] = *da * dx[i__ -1];
		/* L30: */
    }
    if (*PARAM_n < 5) {
		return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *PARAM_n;
    for (i__ = mp1; i__ <= i__1; i__ += 5) {
		dx[i__ -1] = *da * dx[i__ -1];
		dx[i__ + 1 -1] = *da * dx[i__ + 1 -1];
		dx[i__ + 2 -1] = *da * dx[i__ + 2 -1];
		dx[i__ + 3 -1] = *da * dx[i__ + 3 -1];
		dx[i__ + 4 -1] = *da * dx[i__ + 4 -1];
		/* L50: */
    }
    return 0;
} /* dscal_ */

/* DECK DDOT */
__device__ double  ddot_(int * __restrict__ PARAM_n, double  * __restrict__dx, int incx, double  * __restrict__ dy, int incy, struct cuLsodaCommonBlock * __restrict__ common)
{
    /* System generated locals */
    int i__1 = 0;
	int i__2 = 0;
    double  ret_val = 0.;
	
    /* Local variables */
     int i__ = 0;
	int m = 0;
	int ix = 0;
	int iy = 0;
	int ns = 0;
	int mp1 = 0;
	
	/* ***BEGIN PROLOGUE  DDOT */
	/* ***PURPOSE  Compute the inner product of two vectors. */
	/* ***CATEGORY  D1A4 */
	/* ***TYPE      DOUBLE PRECISION (SDOT-S, DDOT-D, CDOTU-C) */
	/* ***KEYWORDS  BLAS, INNER PRODUCT, LINEAR ALGEBRA, VECTOR */
	/* ***AUTHOR  Lawson, C. L., (JPL) */
	/*           Hanson, R. J., (SNLA) */
	/*           Kincaid, D. R., (U. of Texas) */
	/*           Krogh, F. T., (JPL) */
	/* ***DESCRIPTION */
	
	/*                B L A S  Subprogram */
	/*    Description of Parameters */
	
	/*     --Input-- */
	/*        N  number of elements in input vector(s) */
	/*       DX  double  precision vector with N elements */
	/*     INCX  storage spacing between elements of DX */
	/*       DY  double  precision vector with N elements */
	/*     INCY  storage spacing between elements of DY */
	
	/*     --Output-- */
	/*     DDOT  double  precision dot product (zero if N .LE. 0) */
	
	/*     Returns the dot product of double  precision DX and DY. */
	/*     DDOT = sum for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY), */
	/*     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is */
	/*     defined in a similar way using INCY. */
	
	/* ***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. */
	/*                 Krogh, Basic linear algebra subprograms for Fortran */
	/*                 usage, Algorithm No. 539, Transactions on Mathematical */
	/*                 Software 5, 3 (September 1979), pp. 308-323. */
	/* ***ROUTINES CALLED  (NONE) */
	/* ***REVISION HISTORY  (YYMMDD) */
	/*   791001  DATE WRITTEN */
	/*   890831  Modified array declarations.  (WRB) */
	/*   890831  REVISION DATE from Version 3.2 */
	/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
	/*   920310  Corrected definition of LX in DESCRIPTION.  (WRB) */
	/*   920501  Reformatted the REFERENCES section.  (WRB) */
	/* ***END PROLOGUE  DDOT */
	/* ***FIRST EXECUTABLE STATEMENT  DDOT */
   
	 /* Parameter adjustments */
    //--dy;
    //--dx;
	
    /* Function Body */
    ret_val = 0.;
    if (*PARAM_n <= 0) {
		return ret_val;
    }
    if (incx == incy) {
		if ((i__1 = incx - 1) < 0) {
			goto L5;
		} else if (i__1 == 0) {
			goto L20;
		} else {
			goto L60;
		}
    }
	
	/*     Code for unequal or nonpositive increments. */
	
L5:
    ix = 1;
    iy = 1;
    if (incx < 0) {
		ix = (-(*PARAM_n) + 1) * incx + 1;
    }
    if (incy < 0) {
		iy = (-(*PARAM_n) + 1) * incy + 1;
    }
    i__1 = *PARAM_n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		ret_val += dx[ix -1] * dy[iy -1];
		ix += incx;
		iy += incy;
		/* L10: */
    }
    return ret_val;
	
	/*     Code for both increments equal to 1. */
	
	/*     Clean-up loop so remaining vector length is a multiple of 5. */
	
L20:
    m = *PARAM_n % 5;
    if (m == 0) {
		goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
		ret_val += dx[i__ -1] * dy[i__ -1];
		/* L30: */
    }
    if (*PARAM_n < 5) {
		return ret_val;
    }
L40:
    mp1 = m + 1;
    i__1 = *PARAM_n;
    for (i__ = mp1; i__ <= i__1; i__ += 5) 
	{
		ret_val = ret_val + dx[i__ -1] * dy[i__ -1] + dx[i__ + 1 -1] * dy[i__ + 1 -1] + 
		dx[i__ + 2 -1] * dy[i__ + 2 -1] + dx[i__ + 3 -1] * dy[i__ + 3 -1] + dx[i__ + 4 -1] * dy[i__ + 4 -1];
		/* L50: */
    }
    return ret_val;
	
	/*     Code for equal, positive, non-unit increments. */
	
L60:
    ns = *PARAM_n * incx;
    i__1 = ns;
    i__2 = incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) 
	{
		ret_val += dx[i__ -1] * dy[i__ -1];
		/* L70: */
    }
    return ret_val;
} /* ddot_ */

__device__ double  d_sign(double  * __restrict__ a, double  * __restrict__ b)
{
	double  x = 0.;
	x = (*a >= 0 ? *a : - *a);
	return( *b >= 0 ? x : -x);
}

__host__ __device__ void cuLsodaCommonBlockInit(struct cuLsodaCommonBlock * __restrict__ common)
{
	
	/* Common block initialization */
	for (int bugger = 0; bugger < 13; bugger++)
	{
		common->CM_el[bugger] = 0.;
	}
	for (int bugger = 0; bugger < 156; bugger++)
	{
		common->CM_elco[bugger] = 0.;
	}
	for (int bugger = 0; bugger < 36; bugger++)
	{
		common->CM_tesco[bugger] = 0.;
	}
	for (int bugger = 0; bugger < 218; bugger++)
	{
		common->CM_rls[bugger] = 0.;
	}
	for (int bugger = 0; bugger < 12; bugger++)
	{
		common->CM_cm1[bugger] = 0.;
	}
	for (int bugger = 0; bugger < 5; bugger++)
	{
		common->CM_cm2[bugger] = 0.;
	}
	for (int bugger = 0; bugger < 22; bugger++)
	{
		common->CM_rlsa[bugger] = 0.;
	}
	for (int bugger = 0; bugger < 37; bugger++)
	{
		common->CM_ils[bugger] = 0;
	}
		for (int bugger = 0; bugger < 9; bugger++)
	{
		common->CM_ilsa[bugger] = 0;
	}
	double  smThing[12] = { .5,.575,.55,.45,.35,.25,.2,.15,.1,.075,.05,.025 };
	for(int bob = 0; bob <12; bob ++)
	{
		common->CM_sm1[bob] = smThing[bob];
		
	}
	
	// initialize double s in the common block to zero
	common->CM_conit = 0.;
	common->CM_crate = 0.;
	common->CM_ccmax = 0.;
	common->CM_el0 = 0.;
	common->CM_h__ = 0.;
	common->CM_hmin = 0.;
	common->CM_hmxi = 0.;
	common->CM_hu = 0.;
	common->CM_rc = 0.;
	common->CM_tn = 0.;
	common->CM_uround = 0.;
	common->CM_pdest = 0.;
	common->CM_pdlast = 0.;
	common->CM_ratio = 0.;
	common->CM_hold = 0.;
	common->CM_rmax = 0.;
	common->CM_tsw = 0.;
	common->CM_pdnorm = 0.;

	// initialize ints in common block to zero
	common->CM_init = 0;
	common->CM_mxstep = 0;
	common->CM_mxhnil = 0;
	common->CM_nhnil = 0;
	common->CM_nslast = 0;
	common->CM_nyh = 0;
	common->CM_icf = 0;
	common->CM_ierpj = 0;
	common->CM_iersl = 0;
	common->CM_jcur = 0;
	common->CM_jstart = 0;
	common->CM_kflag = 0;
	common->CM_l = 0;
	common->CM_lyh = 0;
	common->CM_lewt = 0;
	common->CM_lacor = 0;
	common->CM_lsavf = 0;
	common->CM_lwm = 0;
	common->CM_liwm = 0;
	common->CM_meth = 0;
	common->CM_miter = 0;
	common->CM_maxord = 0;
	common->CM_maxcor = 0;
	common->CM_msbp = 0;
	common->CM_mxncf = 0;
	common->CM_n = 0;
	common->CM_nq = 0;
	common->CM_nst = 0;
	common->CM_nfe = 0;
	common->CM_nje = 0;
	common->CM_nqu = 0;
	common->CM_ialth = 0;
	common->CM_ipup = 0;
	common->CM_lmax = 0;
	common->CM_nqnyh = 0;
	common->CM_nslp = 0;
	common->CM_insufr = 0;
	common->CM_insufi = 0;
	common->CM_ixpr = 0;
	common->CM_jtyp = 0;
	common->CM_mused = 0;
	common->CM_mxordn = 0;
	common->CM_mxords = 0;
	common->CM_icount = 0;
	common->CM_irflag = 0;
	
	/* End Common Block initialization */

}


template<typename Fex, typename Jex>
__global__ void dlsoda_kernel(Fex fex, 
                              int32_t *     __restrict__ neq, 
                              double  * __restrict__ y, 
                              double  * __restrict__ t, 
                              double  * __restrict__ tout, 
                              int32_t *     __restrict__ itol, 
                              double  * __restrict__ rtol, 
                              double  * __restrict__ atol, 
                              int32_t *     __restrict__ itask, 
                              int32_t *     __restrict__ istate, 
                              int32_t *     __restrict__ iopt, 
                              double  * __restrict__ rwork, 
                              int32_t *     __restrict__ lrw, 
                              int32_t *     __restrict__ iwork, 
                              int32_t *     __restrict__ liw, 
                              Jex jac, 
                              int32_t *     __restrict__ jt, 
                              struct cuLsodaCommonBlock * __restrict__ common, 
                              int32_t *  __restrict__ err, 
                              const int32_t n_threads,
                              const int32_t NEQ)
{
        const int32_t leni = 20+NEQ;
        /* RWORK  = real work array of length at least: */
	/*             22 + NEQ * MAX(16, NEQ + 9). */
        const int32_t lenr = 22+NEQ*max(16,NEQ+9);
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
//	printf("Thread ID: %d\tProbsize: %d\n",me,probSize);
	if(tid < n_threads){
//	printf("neq: %d\ty[0]: %f\ty[1]: %f\ty[2]: %f\ty[3]: %f\tt: %f\ttout: %f\n",neq[me],y[4*me],y[4*me+1],y[4*me+2],y[4*me+3],t[me],tout[me]);
	err[tid] = dlsoda_(fex, 
                           &neq[tid], 
                           &y[NEQ*tid], 
                           &t[tid], 
                           &tout[tid], 
                           &itol[tid], 
                           &rtol[tid], 
                           &atol[tid], 
                           &itask[tid], 
                           &istate[tid], 
                           &iopt[tid], 
                           &rwork[lenr*tid], 
                           &lrw[tid], 
                           &iwork[leni*tid], 
                           &liw[tid], 
                           jac, 
                           &jt[tid], 
                           &common[tid]);
	
	}
	__syncthreads();
}

#include "GMS_cuda_memops.cuh"
#if (PROFILE_HOST_TO_DEVICE) == 1
#include <immintrin.h> //rdtscp
#endif

static const uint64_t rdtscp_cost = 42; // Skylake uarch

template<typename Fex, typename Jex>
void dlsoda_single_gpu(Fex fex,
                int32_t  * __restrict__ neq,
                double  *     __restrict__ y,
                double  *     __restrict__ t, 
                double  *     __restrict__ tout, 
                int32_t *  __restrict__ itol, 
                double  *     __restrict__ rtol, 
                double  *     __restrict__ atol, 
                int32_t *  __restrict__ itask, 
                int32_t *  __restrict__ istate, 
                int32_t *  __restrict__ iopt, 
                double  *     __restrict__ rwork, 
                int32_t *  __restrict__ lrw, 
                int32_t *  __restrict__ iwork, 
                int32_t *  __restrict__ liw, 
                Jex jac, 
                int32_t *  __restrict__ jt, 
                struct cuLsodaCommonBlock * __restrict__ common, 
                int32_t * __restrict__ err, 
                const int32_t n_threads,
                const int32_t NEQ,
                int32_t * ierr,
                cudaError_t * __restrict__ cuerr,
                uint64_t * __restrict tsc_delta ) {
 
       *ierr = 0;
       double  * __restrict__ d_t    = NULL;
       double  * __restrict__ d_y    = NULL;
       int32_t * __restrict__ d_jt   = NULL;
       int32_t * __restrict__ d_neq  = NULL;
       int32_t * __restrict__ d_liw  = NULL;
       int32_t * __restrict__ d_lrw  = NULL;
       double *  __restrict__ d_atol = NULL;
       int32_t * __restrict__ d_itol = NULL;
       int32_t * __restrict__ d_iopt = NULL;
       double  * __restrict__ d_rtol = NULL;
       int32_t * __restrict__ d_iout = NULL;
       double  * __restrict__ d_tout = NULL;
       int32_t * __restrict__ d_itask= NULL;
       int32_t * __restrict__ d_iwork= NULL;
       double  * __restrict__ d_rwork= NULL;
       int32_t * __restrict__ d_istate=NULL;
       struct cuLsodaCommonBlock * __restrict__ d_common = NULL;
       int32_t * __restrict__ d_err  = NULL;
#if (PROFILE_HOST_TO_DEVICE) == 1
       volatile uint64_t dummy1;
       volatile uint32_t dummy2;
       volatile uint64_t tsc_start,tsc_end;
       volatile uint32_t coreid;
#endif
       const int32_t leni = 20+NEQ;
        /* RWORK  = real work array of length at least: */
	/*             22 + NEQ * MAX(16, NEQ + 9). */
       const int32_t lenr = 22+NEQ*max(16,NEQ+9);
       const int32_t threads_blocks = 32;
       const int32_t blocks_grid    = (n_threads+threads_blocks-1)/threads_blocks;
       cudaError_t status;
       int32_t merr = 0;
       alloc_double_gpu(&d_t[0],(size_t)n_threads,&merr);
       alloc_double_gpu(&d_y[0],(size_t)(n_threads*NEQ),&merr);
       alloc_int32_gpu(&d_jt[0],(size_t)n_threads,&merr);
       alloc_int32_gpu(&d_neq[0],(size_t)n_threads,&merr);
       alloc_int32_gpu(&d_liw[0],(size_t)n_threads,&merr);
       alloc_int32_gpu(&d_lrw[0],(size_t)n_threads,&merr);
       alloc_double_gpu(&d_atol[0],(size_t)n_threads,&merr);
       alloc_int32_gpu(&d_itol[0],(size_t)n_threads,&merr);
       alloc_int32_gpu(&d_iopt[0],(size_t)n_threads,&merr);
       alloc_double_gpu(&d_rtol[0],(size_t)n_threads,&merr);
       alloc_int32_gpu(&d_iout[0],(size_t)n_threads,&merr);
       alloc_double_gpu(&d_tout[0],(size_t)n_threads,&merr);
       alloc_int32_gpu(&d_itask[0],(size_t)n_threads,&merr);
       alloc_int32_gpu(&d_iwork[0],(size_t)(n_threads*leni),&merr);
       alloc_double_gpu(&d_rwork[0],(size_t)(n_threads*lenr),&merr);
       alloc_int32_gpu(&d_istate[0],(size_t)n_threads,&merr);
       GMS_CUDA_DEBUG_CHECK(cudaMalloc((void**)&d_common,sizeof(struct cuLsodaCommonBlock)*n_threads));
       alloc_int32_gpu(&d_err[0],(size_t)n_threads,&merr);
       copy_double_cpu_to_gpu(d_t,t,(size_t)n_threads,&merr);
       copy_double_cpu_to_gpu(d_y,y,(size_t)(n_threads*NEQ),&merr);
       copy_int32_cpu_to_gpu(d_jt,jt,(size_t)n_threads,&merr);
       copy_int32_cpu_to_gpu(d_neq,neq,(size_t)n_threads,&merr);
       copy_int32_cpu_to_gpu(d_liw,liw,(size_t)n_threads,&merr);
       copy_int32_cpu_to_gpu(d_lrw,lrw,(size_t)n_threads,&merr);
       copy_double_cpu_to_gpu(d_atol,atol,(size_t)n_threads,&merr);
       copy_int32_cpu_to_gpu(d_itol,itol,(size_t)n_threads,&merr);
       copy_int32_cpu_to_gpu(d_iopt,iopt,(size_t)n_threads,&merr);
       copy_double_cpu_to_gpu(d_rtol,rtol,(size_t)n_threads,&merr);
       copy_int32_cpu_to_gpu(d_iout,iout,(size_t)n_threads,&merr);
       copy_double_cpu_to_gpu(d_tout,tout,(size_t)n_threads,&merr);
       copy_int32_cpu_to_gpu(d_itask,itask,(size_t)n_threads,&merr);
       copy_int32_cpu_to_gpu(d_iwork,iwork,(size_t)(n_threads*leni),&merr);
       copy_double_cpu_to_gpu(d_rwork,rwork,(size_t)(n_threads*lenr),&merr);
       copy_int32_cpu_to_gpu(d_istate,istate,(size_t)n_threads,&merr);
       GMS_CUDA_DEBUG_CHECK(cudaMemcpy(d_common,common,sizeof(struct cuLsodaCommonBlock)*n_threads, cudaMemcpyHostToDevice));
       copy_int32_cpu_to_gpu(d_err,err,(size_t)n_threads,&merr);
#if (PROFILE_HOST_TO_DEVICE) == 1
      dummy1    = __rdtscp(&dummy2);
      __asm__("lfence");
      tsc_start = __rdtscp(&coreid);
      __asm__("lfence");
#endif
      dlsoda_kernel<<<blocks_grid,threads_blocks>>>(fex,d_neq,d_y,d_t,d_tout,d_itol,
                                                    d_rtol,d_atol,d_itask,d_istate,
                                                    d_iopt,d_rwork,d_lrw,d_iwork,d_liw,
                                                    jex,d_jt,d_common,d_err,n_threads);

#if (PROFILE_HOST_TO_DEVICE) == 1
      __asm__("lfence");
      tsc_end     = __rdtscp(&coreid);
      *tsc_delta  = tsc_end-tsc_start-rdtscp_cost;
      __asm__("lfence");
      printf("dlsoda_kernel executed in: %llu reference cycles\n",*tsc_delta);
#endif    
      GMS_CUDA_DEBUG_CHECK(cudaDeviceSynchronize());
      copy_double_gpu_to_cpu(d_t,t,(size_t)n_threads,&merr); 
      copy_double_gpu_to_cpu(d_y,y,(size_t)(n_threads*NEQ),&merr); 
      copy_int32_gpu_to_cpu(d_jt,jt,(size_t)n_threads,&merr);       
      copy_int32_gpu_to_cpu(d_neq,neq,(size_t)n_threads,&merr);      
      copy_int32_gpu_to_cpu(d_liw,liw,(size_t)n_threads,&merr); 
      copy_int32_gpu_to_cpu(d_lrw,lrw,(size_t)n_threads,&merr); 
      copy_double_gpu_to_cpu(d_atol,atol,(size_t)n_threads,&merr);      
      copy_int32_gpu_to_cpu(d_itol,itol,(size_t)n_threads,&merr); 
      copy_int32_gpu_to_cpu(d_iopt,iopt,(size_t)n_threads,&merr);    
      copy_double_gpu_to_cpu(d_rtol,rtol,(size_t)n_threads,&merr); 
      copy_double_gpu_to_cpu(d_tout,tout,(size_t)n_threads,&merr);  
      copy_int32_gpu_to_cpu(d_itask,itask,(size_t)n_threads,&merr);
      copy_int32_gpu_to_cpu(d_iwork,iwork,(size_t)(n_threads*leni),&merr);
      copy_double_gpu_to_cpu(d_rwork,rwork,(size_t)(n_threads*lenr),&merr);
      copy_int32_gpu_to_cpu(d_istate,istate,(size_t)n_threads,&merr);
      GMS_CUDA_DEBUG_CHECK(cudaMemcpy(d_common,common,sizeof(struct cuLsodaCommonBlock)*n_threads, cudaMemcpyDeviceToHost));
      copy_int32_gpu_to_cpu(d_err,err,(size_t)n_threads,&merr); 
      cudaFree(d_t);
      cudaFree(d_y);   
      cudaFree(d_jt); 
      cudaFree(d_neq);  
      cudaFree(d_liw);
      cudaFree(d_lrw); 
      cudaFree(d_atol); 
      cudaFree(d_itol);
      cudaFree(d_iopt);  
      cudaFree(d_rtol); 
      cudaFree(d_tout); 
      cudaFree(d_itask); 
      cudaFree(d_iwork);  
      cudaFree(d_rwork);  
      cudaFree(d_istate); 
      cudaFree(d_common); 
      cudaFree(d_err); 
      *cuerr = status;
      return;

      Error:
      if(d_t)    cudaFree(d_t);
      if(d_y)    cudaFree(d_y);   
      if(d_jt)   cudaFree(d_jt); 
      if(d_neq)  cudaFree(d_neq);  
      if(d_liw)  cudaFree(d_liw);
      if(d_lrw)  cudaFree(d_lrw); 
      if(d_atol) cudaFree(d_atol); 
      if(d_itol) cudaFree(d_itol);
      if(d_iopt) cudaFree(d_iopt);  
      if(d_rtol) cudaFree(d_rtol); 
      if(d_tout) cudaFree(d_tout); 
      if(d_itask)cudaFree(d_itask); 
      if(d_iwork)cudaFree(d_iwork);  
      if(d_rwork)cudaFree(d_rwork);  
      if(d_istate)cudaFree(d_istate); 
      if(d_common)cudaFree(d_common); 
      if(d_err)cudaFree(d_err); 
      *cuerr = status;
      return;
                             
}

#endif /*__GMS_CULSODA_CUH__*/
