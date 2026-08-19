//#define _USE_MATH_DEFINES
#ifndef __GMS_CQUADPACK_H__
#define __GMS_CQUADPACK_H__ 010220221057


#include <math.h>
#include <float.h>
#include "GMS_config.h"


#define uflow     DBL_MIN
#define oflow     DBL_MAX
#define epmach     DBL_EPSILON
#define LIMIT     500
#define MAXP1     21
#ifdef M_PI
#define Pi      M_PI
#else
#define Pi      3.14159265358979323846
#endif
#define COSINE     1
#define SINE    2

#ifndef FALSE
#define FALSE   0
#endif
#ifndef TRUE
#define TRUE    1
#endif
#ifndef min
#define min(a,b)    (((a) < (b)) ? (a) : (b))
#endif
#ifndef max
#define max(a,b)    (((a) > (b)) ? (a) : (b))
#endif



typedef double(*dq_function_type)(double, void * __restrict);

/* Integration routines */
/* Gauss-Kronrod for integration over finite range. */
 double G_K15(dq_function_type f,double a,double b,double * __restrict  abserr,
    double * __restrict  resabs, double * __restrict  resasc, void * __restrict   user_data)
                                                                                              __ATTR_HOT__
                                                                                              __ATTR_ALIGN__(32);
 double G_K21(dq_function_type f, double a, double b, double * __restrict  abserr,
    double * __restrict  resabs, double * __restrict  resasc, void * __restrict   user_data)  __ATTR_HOT__
                                                                                              __ATTR_ALIGN__(32);
 double G_K31(dq_function_type f, double a, double b, double * __restrict  abserr,
    double * __restrict  resabs, double * __restrict  resasc, void * __restrict   user_data)  __ATTR_HOT__
                                                                                              __ATTR_ALIGN__(32);
 double G_K41(dq_function_type f, double a, double b, double * __restrict  abserr,
    double * __restrict  resabs, double * __restrict  resasc, void * __restrict   user_data)  __ATTR_HOT__
                                                                                              __ATTR_ALIGN__(32);
 double G_K51(dq_function_type f, double a, double b, double * __restrict  abserr,
    double * __restrict  resabs, double * __restrict  resasc, void * __restrict   user_data)  __ATTR_HOT__
                                                                                              __ATTR_ALIGN__(32);
 double G_K61(dq_function_type f, double a, double b, double * __restrict  abserr,
    double * __restrict  resabs, double * __restrict  resasc, void * __restrict   user_data)  __ATTR_HOT__
                                                                                              __ATTR_ALIGN__(32);

/* Gauss-Kronrod for integration over infinite range. */
 double G_K15I(dq_function_type f, double boun, int inf, double a, double b,
    double * __restrict  abserr,double * __restrict  resabs, 
    double * __restrict  resasc, void * __restrict   user_data) __ATTR_HOT__
                                                                __ATTR_ALIGN__(32);

/* Gauss-Kronrod for integration of weighted function. */
/*
 __ATTR_HOT__
 __ATTR_ALIGN__(32)
 double G_K15W(dq_function_type f, double w(), double p1, double p2, double p3,
    double p4,int kp,double a,double b,double * __restrict  abserr,
    double * __restrict  resabs, double * __restrict  resasc, void * __restrict   user_data);
*/
 __ATTR_HOT__
 __ATTR_ALIGN__(32)
 double G_K15W(dq_function_type f,double (*w)(double, double, double, double, double, int),
    double p1, double p2, double p3,
    double p4,int kp,double a,double b,double * __restrict  abserr,
    double * __restrict  resabs, double * __restrict  resasc, void * __restrict   user_data);

 double dqext(int * __restrict  n, double epstab [], double * __restrict  abserr,
    double res3la[],int * __restrict  nres) __ATTR_HOT__
                                            __ATTR_ALIGN__(32);
 void dqsort(int limit, int last, int * __restrict  maxerr, double * __restrict  ermax,
    double elist[],int iord[],int * __restrict  nrmax)__ATTR_HOT__
                                                      __ATTR_ALIGN__(32);
 double dqagi(dq_function_type f, double bound, int inf, double epsabs,
    double epsrel,double * __restrict  abserr,int * __restrict  neval,int * __restrict  ier, void * __restrict   user_data) __ATTR_HOT__
                                                                                                                            __ATTR_ALIGN__(32);
 double dqags(dq_function_type f, double a, double b, double epsabs,
    double epsrel,double * __restrict  abserr,int * __restrict  neval,int * __restrict  ier, void * __restrict   user_data) __ATTR_HOT__
                                                                                                                            __ATTR_ALIGN__(32);
 double dqagp(dq_function_type f, double a, double b, int npts2, double * __restrict  points,
    double epsabs,double epsrel,double * __restrict  abserr,int * __restrict  neval,int * __restrict  ier, void * __restrict   user_data) __ATTR_HOT__
                                                                                                                                          __ATTR_ALIGN__(32);
 double dqng(dq_function_type f, double a, double b, double epsabs, double epsrel,
    double * __restrict  abserr,int * __restrict  neval,int * __restrict  ier, void * __restrict   user_data) __ATTR_HOT__
                                                                                                              __ATTR_ALIGN__(32);
 double dqag(dq_function_type f, double a, double b, double epsabs, double epsrel,
    int irule,double * __restrict  abserr,int * __restrict  neval,int * __restrict  ier, void * __restrict   user_data) __ATTR_HOT__
                                                                                                                        __ATTR_ALIGN__(32);
 double dqage(dq_function_type f, double a, double b, double epsabs, double epsrel,
    int irule,double * __restrict  abserr,int * __restrict  neval,int * __restrict  ier,int * __restrict  last, void * __restrict   user_data) __ATTR_HOT__
                                                                                                                                               __ATTR_ALIGN__(32);
 double dqwgtc(double x, double c, double p2, double p3, double p4,
    int kp) __ATTR_HOT__
            __ATTR_ALIGN__(32);
 double dqwgto(double x, double omega, double p2, double p3, double p4,
    int integr) __ATTR_HOT__
                __ATTR_ALIGN__(32);
 double dqwgts(double x, double a, double b, double alpha, double beta,
    int integr) __ATTR_HOT__
                __ATTR_ALIGN__(32);
 void dqcheb(double * __restrict  x, double * __restrict  fval, double * __restrict  cheb12, double * __restrict  cheb24) __ATTR_HOT__
                                                                                                                          __ATTR_ALIGN__(32);
 double dqc25o(dq_function_type f, double a, double b, double omega, int integr,
    int nrmom,int maxp1,int ksave,double * __restrict  abserr,int * __restrict  neval,
    double * __restrict  resabs,double * __restrict  resasc,int * __restrict  momcom,double ** chebmo, void * __restrict   user_data) __ATTR_HOT__
                                                                                                                                                  __ATTR_ALIGN__(32);
 double dqfour(dq_function_type f, double a, double b, double omega, int integr,
    double epsabs,double epsrel,int icall,int maxp1,
    double * __restrict  abserr,int * __restrict  neval,int * __restrict  ier,int * __restrict  momcom,
    double **  chebmo, void * __restrict   user_data) __ATTR_HOT__
                                                                __ATTR_ALIGN__(32);
 double dqawfe(dq_function_type f, double a, double omega, int integr, double epsabs,
    int limlst,int maxp1,double * __restrict  abserr,int * __restrict  neval,int * __restrict  ier,
    double * __restrict  rslst,double * __restrict  erlist,int * __restrict  ierlst,double ** chebmo, void * __restrict   user_data) __ATTR_HOT__
                                                                                                                                                 __ATTR_ALIGN__(32);
 double dqawf(dq_function_type f, double a, double omega, int integr, double epsabs,
    double * __restrict  abserr,int * __restrict  neval,int * __restrict  ier, void * __restrict   user_data) __ATTR_HOT__
                                                                                                              __ATTR_ALIGN__(32);
 double dqawo(dq_function_type f, double a, double b, double omega, int integr, double epsabs,
    double epsrel,double * __restrict  abserr,int * __restrict  neval,int * __restrict  ier, void * __restrict   user_data) __ATTR_HOT__
                                                                                                                            __ATTR_ALIGN__(32);
 double dqaws(dq_function_type f, double a, double b, double alfa, double beta, int wgtfunc,
    double epsabs,double epsrel,double * __restrict  abserr,int * __restrict  neval,int * __restrict  ier, void * __restrict   user_data) __ATTR_HOT__
                                                                                                                                          __ATTR_ALIGN__(32);
 double dqawse(dq_function_type f, double a, double b, double alfa, double beta,
    int wgtfunc,double epsabs,double epsrel,double * __restrict  abserr,
    int * __restrict  neval,int * __restrict  ier, void * __restrict   user_data) __ATTR_HOT__
                                                                                  __ATTR_ALIGN__(32);
 void dqmomo(double alfa, double beta, double ri [], double rj [], double rg [],
    double rh[],int wgtfunc) __ATTR_HOT__
                             __ATTR_ALIGN__(32);
 double dqc25s(dq_function_type f, double a, double b, double bl, double br, double alfa,
    double beta,double ri[],double rj[],double rg[],double rh[],
    double * __restrict  abserr,double * __restrict  resasc,int wgtfunc,int * __restrict  nev, void * __restrict   user_data) __ATTR_HOT__
                                                                                                                              __ATTR_ALIGN__(32);
 double dqc25c(dq_function_type f, double a, double b, double c, double * __restrict  abserr,
    int * __restrict  krul,int * __restrict  neval, void * __restrict   user_data) __ATTR_HOT__
                                                                                   __ATTR_ALIGN__(32);
 double dqawc(dq_function_type f, double a, double b, double c, double epsabs,
    double epsrel,double * __restrict  abserr,int * __restrict  neval,int * __restrict  ier, void * __restrict   user_data) __ATTR_HOT__
                                                                                                                            __ATTR_ALIGN__(32);
 double dqawce(dq_function_type f, double a, double b, double c, double epsabs,
    double epsrel,double * __restrict  abserr,int * __restrict  neval,int * __restrict  ier, void * __restrict   user_data) __ATTR_HOT__
                                                                                                                            __ATTR_ALIGN__(32);

 double G_B15(dq_function_type f, double a, double b, double * __restrict  abserr,
    double * __restrict  resabs, double * __restrict  resasc, void * __restrict   user_data) __ATTR_HOT__
                                                                                             __ATTR_ALIGN__(32);

#endif /*__GMS_CQUADPACK_H__*/
