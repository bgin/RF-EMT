

#include "GMS_cquadpack.h"


/* dqag.c -- modified version of QUADPACK routine DQAG.
 * (C)1999, C. Bond. All right reserved.
 *
 * There are no changes to the basic computational method. Only
 * the temporary storage strategy is changed to utilize the
 * local stack at the appropriate level. This reduces the
 * need for memory allocation of arrays at higher levels and
 * the resulting passing of memory pointers down the line.
 *
 */


/* DQAG - Approximation to definite integral. (From QUADPACK)
 *
 *  Calls DQAGE with appropriate parameters assigned.
 *
 * PARAMETERS:
 *
 *    f() - double precision function to be integrated.
 *
 *    a - lower limit of integration.
 *
 *    b - upper limit of integration.
 *
 *    epsabs - absolute accuracy requested.
 *
 *    epsrel - relative accuracy requested.
 *
 *    irule - integration rule to be used as follows:
 *        irule = 1 -- G_K 7-15
 *        irule = 2 -- G_K 10-21
 *        irule = 3 -- G_K 15-31
 *        irule = 4 -- G_K 20-41
 *        irule = 5 -- G_K 25-51
 *        irule = 6 -- G_K 30-61
 */
double dqag(dq_function_type f,double a,double b,double epsabs,
    double epsrel,int irule,double *abserr,int *neval,int *ier, void* user_data)
{
    double result;
    int last;

    result = dqage(f,a,b,epsabs,epsrel,irule,abserr,neval,ier,&last, user_data);

    return result;
}



/* DQAGE - Approximation to definite integral. (From QUADPACK)
 *
 *    Allows user's choice of Gauss-Kronrod integration rule.
 *
 * PARAMETERS:
 *
 *    f() - double precision function to be integrated.
 *
 *    a - lower limit of integration.
 *
 *    b - upper limit of integration.
 *
 *    epsabs - absolute accuracy requested.
 *
 *    epsrel - relative accuracy requested.
 *
 *    irule - integration rule to be used as follows:
 *        irule = 1 -- G_K 7-15
 *        irule = 2 -- G_K 10-21
 *        irule = 3 -- G_K 15-31
 *        irule = 4 -- G_K 20-41
 *        irule = 5 -- G_K 25-51
 *        irule = 6 -- G_K 30-61
 *
 *    limit - maximum number of subintervals.
 */
double dqage(dq_function_type f,double a,double b,double epsabs,double epsrel,
    int irule,double *abserr,int *neval,int *ier,int *last, void* user_data)
{
    double area,area1,area2,area12,a1,a2,b1,b2,c,defabs;
    double defab1,defab2,errbnd,errmax,error1,error2;
    double erro12,errsum,resabs,result;
    double alist[LIMIT],blist[LIMIT],rlist[LIMIT],elist[LIMIT];
    int iroff1,iroff2,k,keyf,maxerr,nrmax,iord[LIMIT],limit;

    limit = LIMIT - 1;
    *ier = 0;
    *neval = 0;
    *last = 0;
    result = 0.0;
    *abserr = 0.0;
    alist[0] = a;
    blist[0] = b;
    rlist[0] = 0.0;
    elist[0] = 0.0;
    iord[0] = 0;
    defabs = 0.0;
    resabs = 0.0;
    if ((epsabs < 0.0) && (epsrel < 0.0))
        *ier = 6;
    if (*ier == 6) return result;

/* First approximation to the integral. */
    keyf = irule;
    if (irule <= 0) keyf = 1;
    if (irule >= 7) keyf = 6;
    c = keyf;
    *neval = 0;
    switch (keyf) {
        case 1:
            result = G_K15(f,a,b,abserr,&defabs,&resabs, user_data);
            break;
        case 2:
            result = G_K21(f,a,b,abserr,&defabs,&resabs, user_data);
            break;
        case 3:
            result = G_K31(f,a,b,abserr,&defabs,&resabs, user_data);
            break;
        case 4:
            result = G_K41(f,a,b,abserr,&defabs,&resabs, user_data);
            break;
        case 5:
            result = G_K51(f,a,b,abserr,&defabs,&resabs, user_data);
            break;
        case 6:
            result = G_K61(f,a,b,abserr,&defabs,&resabs, user_data);
            break;
    }
    *last = 0;
    rlist[0] = result;
    elist[0] = *abserr;
    iord[0] = 0;

/* Test on accuracy. */
    errbnd = max(epsabs,epsrel * fabs(result));
    if ((*abserr <= 50.0 * epmach * defabs) && (*abserr > errbnd))
        *ier = 2;
    if (limit == 0) *ier = 1;
    if ((*ier != 0) || ((*abserr <= errbnd) && (*abserr != resabs)) ||
        (*abserr == 0.0)) goto _60;

/* Initialization. */
    errmax = *abserr;
    maxerr = 0;
    area = result;
    errsum = *abserr;
    nrmax = 0;
    iroff1 = 0;
    iroff2 = 0;

/* Main Loop. */
    for (*last = 1; *last <= limit; (*last)++) {
/* Bisect the subinterval with the largest error estimate. */
        a1 = alist[maxerr];
        b1 = 0.5 * (alist[maxerr] + blist[maxerr]);
        a2 = b1;
        b2 = blist[maxerr];
        switch (keyf) {
            case 1:
                area1 = G_K15(f,a1,b1,&error1,&resabs,&defab1, user_data);
                area2 = G_K15(f,a2,b2,&error2,&resabs,&defab2, user_data);
                break;
            case 2:
                area1 = G_K21(f,a1,b1,&error1,&resabs,&defab1, user_data);
                area2 = G_K21(f,a2,b2,&error2,&resabs,&defab2, user_data);
                break;
            case 3:
                area1 = G_K31(f,a1,b1,&error1,&resabs,&defab1, user_data);
                area2 = G_K31(f,a2,b2,&error2,&resabs,&defab2, user_data);
                break;
            case 4:
                area1 = G_K41(f,a1,b1,&error1,&resabs,&defab1, user_data);
                area2 = G_K41(f,a2,b2,&error2,&resabs,&defab2, user_data);
                break;
            case 5:
                area1 = G_K51(f,a1,b1,&error1,&resabs,&defab1, user_data);
                area2 = G_K51(f,a2,b2,&error2,&resabs,&defab2, user_data);
                break;
            case 6:
                area1 = G_K61(f,a1,b1,&error1,&resabs,&defab1, user_data);
                area2 = G_K61(f,a2,b2,&error2,&resabs,&defab2, user_data);
                break;
        }

/* Improve previous approximations to integral and error,
        and test for accuracy. */
        (*neval) += 1;
        area12 = area1 + area2;
        erro12 = error1 + error2;
        errsum = errsum + erro12 - errmax;
        area = area + area12 - rlist[maxerr];
        if ((defab1 != error1) && (defab2 != error2)) {
            if ((fabs(rlist[maxerr]-area12) <= 1.0e-5 * fabs(area12)) &&
                (erro12 >= .99 * errmax))
                    iroff1++;
            if ((*last > 9) && (erro12 > errmax))
                iroff2++;
        }
        rlist[maxerr] = area1;
        rlist[*last] = area2;
        errbnd = max(epsabs,epsrel * fabs(area));
        if (errsum > errbnd)  {

/* Test for roundoff error and eventually set error flag. */
            if ((iroff1 > 6) || (iroff2 > 20))
                *ier = 2;

/* Set error flag in the case that the number of subintervals
    equals the limit. */
            if (*last == limit)
                *ier = 1;

/* Set error flag in the case of bad integrand behavior at a
    point of the integration range. */
            if (max(fabs(a1),fabs(b2)) <= (1.0 + c * 1000.0 * epmach) *
                (fabs(a2)+1.0e4 * uflow))
            *ier = 3;
        }
/* Append the newly-created intervals to the list. */

        if (error2 <= error1) {
            alist[*last] = a2;
            blist[maxerr] = b1;
            blist[*last] = b2;
            elist[maxerr] = error1;
            elist[*last] = error2;
        }
        else {
            alist[maxerr] = a2;
            alist[*last] = a1;
            blist[*last] = b1;
            rlist[maxerr] = area2;
            rlist[*last] = area1;
            elist[maxerr] = error2;
            elist[*last] = error1;
        }

/* Call DQSORT to maintain the descending ordering in the list of
    error estimates and select the subinterval with the
    largest error estimate (to be bisected next). */

        dqsort(limit,*last,&maxerr,&errmax,elist,iord,&nrmax);
        if ((*ier != 0) || (errsum <= errbnd)) break;
    }

/* Compute final result. */

    result = 0.0;
    for (k = 0; k <= *last; k++) {
        result += rlist[k];
    }
    *abserr = errsum;
_60:
    if (keyf != 1)
        *neval = (10 * keyf + 1) * (2 * (*neval) + 1);
    else
        *neval = 30 * (*neval) + 15;

    return result;
}




/* DQAGI - Integration over (semi-) infinite intervals. (From QUADPACK)
 *
 *    Adaptive integration routine which handles functions
 *    to be integrated between -infinity to +infinity, or
 *    between either of those limits and some finite,
 *    real boundary.
 *
 *    The adaptive strategy compares results of integration
 *    over the interval with the sum of results obtained from
 *    integration of bisected interval. Since error estimates
 *    are available from each regional integration, the interval
 *    with the largest error is bisected and new results are
 *    computed. This bisection process is continued until the
 *    error is less than the prescribed limit or convergence
 *    failure is determined.
 *
 *    Note that bisection, in the sense used above, refers to
 *    bisection of the transformed interval.
 *
 * PARAMETERS:
 *
 *    f() - double precision function to be integrated.
 *
 *    bound - optional finite bound on integral.
 *
 *    inf - specifies range of integration as follows:
 *        inf = -1 -- range is from -infinity to bound,
 *        inf =  1 -- range is from bound to +infinity,
 *        inf =  2 -- range is from -infinity to +infinity,
 *                (bound is immaterial in this case).
 *
 *    epsabs - absolute accuracy requested.
 *
 *    epsrel - relative accuracy requested.
 */
double dqagi(dq_function_type f,double bound,int inf,double epsabs,
    double epsrel,double *abserr,int *neval,int *ier, void* user_data)
{
    double abseps,area,area1,area12,area2,a1,a2,b1,b2;
    double boun,correc,defabs,defab1,defab2,dres,erlarg;
    double erlast,errbnd,errmax,error1,error2,erro12;
    double errsum,ertest,resabs,reseps,result,res3la[3];
    double alist[LIMIT],blist[LIMIT],elist[LIMIT],rlist[LIMIT];
    double rlist2[52],small = 0; /* small will be initialized in _80 */

    int id,ierro,iord[LIMIT],iroff1,iroff2,iroff3,jupbnd,k,ksgn;
    int ktmin,last,maxerr,nres,nrmax,numrl2;
    int limit,extrap,noext;

    limit = LIMIT - 1;
/* Test validity of parameters. */
    *ier = 0;
    *neval = 0;
    last = 0;
    result = 0.0;
    *abserr = 0.0;
    alist[0] = 0.0;
    blist[0] = 1.0;
    rlist[0] = 0.0;
    elist[0] = 0.0;
    iord[0] = 0;
    if ((epsabs < 0.0) && (epsrel < 0.0)) *ier = 6;
    if ((inf != 1) && (inf != -1) && (inf != 2)) *ier = 6;
    if (*ier == 6) return result;

/* First approximation to the integral. */
    boun = bound;
    if (inf == 2) boun = 0.0;

    result = G_K15I(f,boun,inf,0.0,1.0,abserr,&defabs,&resabs, user_data);

/* Test on accuracy. */
    last = 0;
    rlist[0] = result;
    elist[0] = *abserr;
    iord[0] = 0;
    dres = fabs(result);
    errbnd = max(epsabs,epsrel*dres);
    if ((*abserr <= 100.0 * epmach * defabs) && (*abserr > errbnd))
        *ier = 2;
    if (limit == 0) *ier = 1;
    if ((*ier != 0) || ((*abserr <= errbnd) && (*abserr != resabs)) ||
        (*abserr == 0.0)) goto _130;

/* Initialization for main loop. */
    rlist2[0] = result;
    errmax = *abserr;
    maxerr = 0;             /* maxerr = 1 */
    area = result;
    errsum = *abserr;
    *abserr = oflow;
    nrmax = 0;
    nres = 0;          /* nres = 0 */
    ktmin = 0;
    numrl2 = 1;            /* numrl2 = 2 */
    extrap = FALSE;
    noext = FALSE;
    ierro = 0;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ksgn = -1;
    if (dres > (1.0 - 50.0 * epmach) * defabs)
        ksgn = 1;

/* Main loop. */
    for (last = 1; last <= limit; last++) {
        a1 = alist[maxerr];
        b1 = 0.5 * (alist[maxerr] + blist[maxerr]);
        a2 = b1;
        b2 = blist[maxerr];
        erlast = errmax;
        area1 = G_K15I(f,boun,inf,a1,b1,&error1,&resabs,&defab1, user_data);
        area2 = G_K15I(f,boun,inf,a2,b2,&error2,&resabs,&defab2, user_data);

/* Improve previous approxminations to integral and error
      and test for accuracy. */
        area12 = area1 + area2;
        erro12 = error1 + error2;
        errsum = errsum + erro12 - errmax;
        area = area + area12 - rlist[maxerr];
        if ((defab1 == error1) || (defab2 == error2)) goto _15;
        if ((fabs(rlist[maxerr] - area12) > 1.0e-5 * fabs(area12))
            || (erro12 < .99 * errmax)) goto _10;
        if (extrap) iroff2++;
        else iroff1++;
_10:
        if ((last > 9) && (erro12 > errmax))    /* last > 10 */
            iroff3++;
_15:
        rlist[maxerr] = area1;
        rlist[last] = area2;
        errbnd = max(epsabs,epsrel * fabs(area));

/* Test for roundoff error and eventually set error flag. */
        if (((iroff1 + iroff2) >= 10) || (iroff3 >= 20))
            *ier = 2;
        if (iroff2 > 5)
            *ier = 3;

/* Set error flag in the case that the number of subintervals
    equals limit. */
        if (last == limit)    /* last == limit */
            *ier = 1;

/* Set error flag in the case of bad integrand behavior at some
    points in the integration range. */
        if (max(fabs(a1),fabs(b2)) <= (1.0 +1000.0 * epmach) *
            (fabs(a2) + 1000.0*uflow))
            *ier = 4;

/* Append the newly-created intervals to the list. */
        if (error2 <= error1) {
            alist[last] = a2;
            blist[maxerr] = b1;
            blist[last] = b2;
            elist[maxerr] = error1;
            elist[last] = error2;
        }
        else {
            alist[maxerr] = a2;
            alist[last] = a1;
            blist[last] = b1;
            rlist[maxerr] = area2;
            rlist[last] = area1;
            elist[maxerr] = error2;
            elist[last] = error1;
       }
/* Call dqsort to maintain the descending ordering in the list of error
    estimates and select the subinterval with nrmax-th largest
    error estimate (to be bisected next). */

        dqsort(limit,last,&maxerr,&errmax,elist,iord,&nrmax);
        if (errsum <= errbnd) goto _115;
        if (*ier != 0) goto _100;
        if (last == 1) goto _80;    /* last == 2 */
        if (noext) continue;  //goto _90;
        erlarg -= erlast;
        if (fabs(b1-a1) > small)
            erlarg += erro12;
        if (extrap) goto _40;

/* Test whether the interval to be bisected next is the smallest interval. */
        if ((fabs(blist[maxerr] - alist[maxerr])) > small)
            goto _90;
        extrap = TRUE;
        nrmax = 1;        /* nrmax = 2 */
_40:
        if ((ierro == 3) || (erlarg <= ertest)) goto _60;

/* The smallest interval has the largest error. Before bisecting, decrease
    the sum of the erors over the larger intervals (erlarg) and
        perform extrapolation.) */
        id = nrmax;
        jupbnd = last;
        if (last > (2 + limit/2))
            jupbnd = limit + 3 - last;
        for (k = id;k <= jupbnd; k++) {
            maxerr = iord[nrmax];
            errmax = elist[maxerr];
            if (fabs(blist[maxerr] - alist[maxerr]) > small)
                goto _90;
            nrmax++;
        }

/* Perform extrapolation. */
_60:
        numrl2++;
        rlist2[numrl2] = area;
        reseps=dqext(&numrl2,rlist2,&abseps,res3la,&nres);
        ktmin++;
        if ((ktmin > 5) && (*abserr < 1.0e-3 * errsum)) *ier = 5;
        if (abseps >= *abserr) goto _70;
        ktmin = 0;
        *abserr = abseps;
        result = reseps;
        correc = erlarg;
        ertest = max(epsabs,epsrel * fabs(reseps));
        if (*abserr <= ertest) goto _100;

/* Prepare bisection of the smallest interval. */
_70:
        if (numrl2 == 0) noext = TRUE;
        if (*ier == 5) goto _100;
        maxerr = iord[0];
        errmax = elist[maxerr];
        nrmax = 0;
        extrap = FALSE;
        small = small * 0.5;
        erlarg = errsum;
        continue;
_80:
        small = .375;
        erlarg = errsum;
        ertest = errbnd;
        rlist2[1] = area;
_90:
        ;
    }                    /* 90: */
_100:
    if (*abserr == oflow) goto _115;
    if ((*ier + ierro) == 0) goto _110;
    if (ierro == 3) *abserr += correc;
    if (*ier == 0) *ier = 3;
    if ((result != 0.0) && (area != 0.0)) goto _105;
    if (*abserr > errsum) goto _115;
    if (area == 0.0) goto _130;
    goto _110;
_105:
    if (*abserr/fabs(result) > errsum/fabs(area)) goto _115;

/* Test on divergence. */
_110:
    if ((ksgn == -1) && (max(fabs(result),fabs(area)) <= defabs * .01))
        goto _130;
    if ((0.01 > result/area) || (result/area > 100.0) ||
        (errsum > fabs(area))) *ier = 6;
    goto _130;

/* Compute global integral. */
_115:
    result = 0.0;
    for (k = 0; k <= last; k++)
        result += rlist[k];
    *abserr = errsum;
_130:
    *neval = 30 * last + 15;
    if (inf == 2) *neval *= 2;
    if (*ier > 2) (*ier)--;
    return result;
}



/* DQAGP - Integration over finite intervals. (From QUADPACK)
 *       Accepts a list of known singularities.
 *
 *    Adaptive integration routine which handles functions
 *    to be integrated between two finite bounds.
 *
 *    The adaptive strategy compares results of integration
 *    over the given interval with the sum of results obtained
 *    from integration over a bisected interval. Since error
 *    estimates are available from each regional integration, the
 *    region with the largest error is bisected and new results
 *    are computed. This bisection process is continued until the
 *    error is less than the prescribed limit or convergence
 *    failure is determined.
 *
 * PARAMETERS:
 *
 *    f() - double precision function to be integrated.
 *
 *    a - lower limit of integration.
 *
 *    b - upper limit of integration.
 *
 *    npts2 - number equal to 2 more than the number of sinularities.
 *
 *    points - vector of dimension npts2, the first (npts2-2) elements
 *         of which are the user provided interior break points.
 *
 *    epsabs - absolute accuracy requested.
 *
 *    epsrel - relative accuracy requested.
 */
double dqagp(dq_function_type f,double a,double b,int npts2,double *points,
    double epsabs,double epsrel,double *abserr,int *neval,int *ier, void* user_data)
{
    double abseps,alist[LIMIT],area,area1,area12,area2;
    double a1,a2,blist[LIMIT],b1,b2,correc,defabs,defab1;
    double defab2,dres,elist[LIMIT],erlarg,erlast,errbnd;
    double errmax,error1,error2,erro12,errsum,ertest,ndin[40];
    double pts[40],resa,resabs,reseps,result,res3la[3];
    double rlist[LIMIT],rlist2[52],sign,temp;

    int i,id,ierro,ind1,ind2,ip1,iord[LIMIT],iroff1,iroff2,iroff3;
    int j,jlow,jupbnd,k,ksgn,ktmin,last,levcur,level[LIMIT],levmax;
    int maxerr,nint,nintp1,npts,nres,nrmax,numrl2,limit,extrap,noext;

    limit = LIMIT - 1;

/* Test validity of parameters. */
    *ier = 0;
    *neval = 0;
    last = 0;
    result = 0.0;
    *abserr = 0.0;
    alist[0] = a;
    blist[0] = b;
    rlist[0] = 0.0;
    elist[0] = 0.0;
    iord[0] = 0;
    level[0] = 0;
    npts = npts2-2;
    if ((npts2 < 2) || (limit < npts) || ((epsabs < 0.0) &&
        (epsrel < 0.0))) *ier = 6;
    if (*ier == 6) goto _999;

/* If any break points are provided, sort them into an ascending sequence. */
    sign = (a < b) ? 1.0 : -1.0;
    pts[0] = min(a,b);
    if (npts == 0) goto _15;
    for (i = 0; i < npts; i++)
        pts[i+1] = points[i];
_15:
    pts[npts+1] = max(a,b);
    nint = npts + 1;
    a1 = pts[0];
    if (npts == 0) goto _40;
    nintp1 = nint + 1;
    for (i = 0; i < nint; i++) {
        ip1 = i + 1;
        for (j = ip1; j < nintp1; j++) {
            if (pts[i] <= pts[j])
                goto _20;
            temp = pts[i];
            pts[i] = pts[j];
            pts[j] = temp;
_20:
            ;
        }
    }
    if ((pts[0] != min(a,b)) || (pts[nintp1-1] != max(a,b)))
        *ier = 6;
    if (*ier == 6)
        goto _999;

/* Compute first integral and error approximations. */
_40:
    resabs = 0.0;
    for (i = 0; i < nint; i++) {
        b1 = pts[i+1];
        area1 = G_K21(f,a1,b1,&error1,&defabs,&resa, user_data);
        *abserr = *abserr + error1;
        result = result + area1;
        ndin[i] = 0;
        if ((error1 == resa) && (error1 != 0.0))
            ndin[i] = 1;
        resabs += defabs;
        level[i] = 0;
        elist[i] = error1;
        alist[i] = a1;
        blist[i] = b1;
        rlist[i] = area1;
        iord[i] = i;
        a1 = b1;
    }
    errsum = 0.0;
    for (i = 0; i < nint; i++) {
        if (ndin[i] == 1)
            elist[i] = *abserr;
        errsum += elist[i];
    }

/* Test on accuracy. */
/*      last = nint; */
    *neval = 21 * nint;
    dres = fabs(result);
    errbnd = max(epsabs,epsrel*dres);
    if ((*abserr <= 100.0 * epmach * resabs) && (*abserr > errbnd))
        *ier = 2;
    if (nint == 0)
        goto _80;
    for (i = 0; i < npts; i++) {
        jlow = i + 1;
        ind1 = iord[i];
        for (j = jlow; j < nint; j++) {
            ind2 = iord[j];
            if (elist[ind1] > elist[ind2])
                goto _60; /* use continue after debugging */
            ind1 = ind2;
            k = j;
_60:
            ;
        }
        if (ind1 == iord[i])
            goto _70;
        iord[k] = iord[i];
        iord[i] = ind1;
_70:
        ;
    }
    if (limit < npts2)
        *ier = 1;
_80:
    if ((*ier != 0) || (*abserr <= errbnd))
        goto _999;

/* Initialization. */
    res3la[0] = 0.0;
    res3la[1] = 0.0;
    res3la[2] = 0.0;
    rlist2[0] = result;
    maxerr = iord[0];
    errmax = elist[maxerr];
    area = result;
    nrmax = 0;
    nrmax = 0;
    nres = -1;            /* nres = 0 */
    numrl2 = 0;            /* numrl2 = 1 */
    ktmin = 0;
    extrap = FALSE;
    noext = FALSE;
    erlarg = errsum;
    ertest = errbnd;
    levmax = 1;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ierro = 0;
    *abserr = oflow;
    ksgn = -1;
    if (dres > (1.0 - 50.0 * epmach) * resabs)
        ksgn = 1;

/* Main loop. */
    for (last = npts2; last <= limit; last++) {

/* Bisect the interval with the nrmax-th largest error estimate. */
        levcur = level[maxerr] + 1;
        a1 = alist[maxerr];
        b1 = 0.5 * (alist[maxerr] + blist[maxerr]);
        a2 = b1;
        b2 = blist[maxerr];
        erlast = errmax;
        area1 = G_K21(f,a1,b1,&error1,&resa,&defab1, user_data);
        area2 = G_K21(f,a2,b2,&error2,&resa,&defab2, user_data);
/* Improve previous approximations to integral and error
      and test for accuracy. */
          *neval += 42;
        area12 = area1 + area2;
        erro12 = error1 + error2;
        errsum = errsum + erro12 - errmax;
        area = area + area12 - rlist[maxerr];
        if ((defab1 == error1) || (defab2 == error2)) goto _95;
        if ((fabs(rlist[maxerr] - area12) > 1.0e-5 * fabs(area12))
            || (erro12 < .99 * errmax)) goto _90;
        if (extrap) iroff2++;
        else iroff1++;
_90:
        if ((last > 9) && (erro12 > errmax))    /* last > 10 */
            iroff3++;
_95:
        level[maxerr] = levcur;
        level[last] = levcur;
        rlist[maxerr] = area1;
        rlist[last] = area2;
        errbnd = max(epsabs,epsrel * fabs(area));

/* Test for roundoff error and eventually set error flag. */
        if (((iroff1 + iroff2) >= 10) || (iroff3 >= 20))
            *ier = 2;
        if (iroff2 > 5)
            *ier = 3;

/* Set error flag in the case that the number of subintervals
    equals limit. */
        if (last == limit)    /* last == limit */
            *ier = 1;

/* Set error flag in the case of bad integrand behavior at some
    points in the integration range. */
        if (max(fabs(a1),fabs(b2)) <= (1.0 +1000.0 * epmach) *
            (fabs(a2) + 1000.0*uflow))
            *ier = 4;

/* Append the newly-created intervals to the list. */
        if (error2 > error1) goto _100;
        alist[last] = a2;
        blist[maxerr] = b1;
        blist[last] = b2;
        elist[maxerr] = error1;
        elist[last] = error2;
        goto _110;
_100:
        alist[maxerr] = a2;
        alist[last] = a1;
        blist[last] = b1;
        rlist[maxerr] = area2;
        rlist[last] = area1;
        elist[maxerr] = error2;
        elist[last] = error1;

/* Call dqsort to maintain the descending ordering in the list of error
    estimates and select the subinterval with nrmax-th largest
    error estimate (to be bisected next). */
_110:
        dqsort(limit,last,&maxerr,&errmax,elist,iord,&nrmax);
        if (errsum <= errbnd) goto _190;
        if (*ier != 0) goto _170;
        if (noext) goto _160;
        erlarg -= erlast;
        if (levcur+1 <= levmax)
            erlarg += erro12;
        if (extrap) goto _120;

/* Test whether the interval to be bisected next is the smallest interval. */
        if ((level[maxerr]+1) <= levmax)
            goto _160;
        extrap = TRUE;
        nrmax = 1;        /* nrmax = 2 */
_120:
        if ((ierro == 3) || (erlarg <= ertest)) goto _140;

/* The smallest interval has the largest error. Before bisecting, decrease
    the sum of the errors over the larger intervals (erlarg) and
        perform extrapolation.) */
        id = nrmax;
        jupbnd = last;
        if (last > (2 + limit/2))
            jupbnd = limit + 3 - last;
        for (k = id;k <= jupbnd; k++) {
            maxerr = iord[nrmax];
            errmax = elist[maxerr];
            if (level[maxerr]+1 <= levmax)
                goto _160;
            nrmax++;
        }

/* Perform extrapolation. */
_140:
        numrl2++;
        rlist2[numrl2] = area;
        if (numrl2 <= 1) goto _155;
        reseps=dqext(&numrl2,rlist2,&abseps,res3la,&nres);
        ktmin++;
        if ((ktmin > 5) && (*abserr < 1.0e-3 * errsum)) *ier = 5;
        if (abseps >= *abserr) goto _150;
        ktmin = 0;
        *abserr = abseps;
        result = reseps;
        correc = erlarg;
        ertest = max(epsabs,epsrel * fabs(reseps));
        if (*abserr <= ertest) goto _170;

/* Prepare bisection of the smallest interval. */
_150:
        if (numrl2 == 0) noext = TRUE;
        if (*ier == 5) goto _170;
_155:
        maxerr = iord[0];
        errmax = elist[maxerr];
        nrmax = 0;
        extrap = FALSE;
        levmax += 1.0;
        erlarg = errsum;
_160:
        ;
    }
_170:
    if (*abserr == oflow) goto _190;
    if ((*ier + ierro) == 0) goto _180;
    if (ierro == 3) *abserr += correc;
    if (*ier == 0) *ier = 3;
    if ((result != 0.0) && (area != 0.0)) goto _175;
    if (*abserr > errsum) goto _190;
    if (area == 0.0) goto _210;
    goto _180;
_175:
    if (*abserr/fabs(result) > errsum/fabs(area)) goto _190;

/* Test on divergence. */
_180:
    if ((ksgn == -1) && (max(fabs(result),fabs(area)) <= defabs * .01))
        goto _210;
    if ((0.01 > result/area) || (result/area > 100.0) ||
        (errsum > fabs(area))) *ier = 6;
    goto _210;

/* Compute global integral. */
_190:
    result = 0.0;
    for (k = 0; k <= last; k++)
        result += rlist[k];
    *abserr = errsum;
_210:
    if (*ier > 2) (*ier)--;
    result = result * sign;
_999:
    return result;
}




/* DQAGS - Integration over finite intervals. (From QUADPACK)
 *
 *    Adaptive integration routine which handles functions
 *    to be integrated between two finite bounds.
 *
 *    The adaptive strategy compares results of integration
 *    over the given interval with the sum of results obtained
 *    from integration over a bisected interval. Since error
 *    estimates are available from each regional integration, the
 *    region with the largest error is bisected and new results
 *    are computed. This bisection process is continued until the
 *    error is less than the prescribed limit or convergence
 *    failure is determined.
 *
 * PARAMETERS:
 *
 *    f() - double precision function to be integrated.
 *
 *    a - lower limit of integration.
 *
 *    b - upper limit of integration.
 *
 *    epsabs - absolute accuracy requested.
 *
 *    epsrel - relative accuracy requested.
 */
double dqags(dq_function_type f,double a,double b,double epsabs,
    double epsrel,double *abserr,int *neval,int *ier, void* user_data)
{
    double abseps,alist[LIMIT],area,area1,area12,area2;
    double a1,a2,blist[LIMIT],b1,b2,correc,defabs,defab1;
    double defab2,dres,elist[LIMIT],erlarg,erlast,errbnd;
    double errmax,error1,error2,erro12,errsum,ertest;
    double resabs,reseps,result,res3la[3],rlist[LIMIT];
    double rlist2[52],small = 0; /* small will be initialized in _80 */

    int id,ierro,iord[LIMIT],iroff1,iroff2,iroff3,jupbnd,k,ksgn;
    int ktmin,last,maxerr,nres,nrmax,numrl2;
    int limit;
    int extrap,noext;

    limit = LIMIT -1;
/* Test validity of parameters. */
    *ier = 0;
    *neval = 0;
    result = 0.0;
    *abserr = 0.0;
    alist[0] = a;
    blist[0] = b;
    rlist[0] = 0.0;
    elist[0] = 0.0;
    if ((epsabs < 0.0) && (epsrel < 0.0)) *ier = 6;
    if (*ier == 6) return result;

/* First approximation to the integral. */
    ierro = 0;
    result = G_K21(f,a,b,abserr,&defabs,&resabs, user_data);

/* Test on accuracy. */
    dres = fabs(result);
    errbnd = max(epsabs,epsrel*dres);
    last = 1;
    rlist[0] = result;
    elist[0] = *abserr;
    iord[0] = 0;
    if ((*abserr <= 100.0 * epmach * defabs) && (*abserr > errbnd))
        *ier = 2;
    if (limit == 0) *ier = 1;
    if ((*ier != 0) || ((*abserr <= errbnd) && (*abserr != resabs)) ||
        (*abserr == 0.0)) goto _140;

/* Initialization. */
    rlist2[0] = result;
    errmax = *abserr;
    maxerr = 0;             /* maxerr = 1 */
    area = result;
    errsum = *abserr;
    *abserr = oflow;
    nrmax = 0;
    nres = 0;          /* nres = 0 */
    numrl2 = 1;            /* numrl2 = 2 */
    ktmin = 0;
    extrap = FALSE;
    noext = FALSE;
    ierro = 0;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ksgn = -1;
    if (dres > (1.0 - 50.0 * epmach) * defabs)
        ksgn = 1;

/* Main loop. */
    for (last = 1; last <= limit; last++) {

/* Bisect the interval with the nrmax-th largest error estimate. */
        a1 = alist[maxerr];
        b1 = 0.5 * (alist[maxerr] + blist[maxerr]);
        a2 = b1;
        b2 = blist[maxerr];
        erlast = errmax;
        area1 = G_K21(f,a1,b1,&error1,&resabs,&defab1, user_data);
        area2 = G_K21(f,a2,b2,&error2,&resabs,&defab2, user_data);

/* Improve previous approximation's to integral and error
      and test for accuracy. */
        area12 = area1 + area2;
        erro12 = error1 + error2;
        errsum = errsum + erro12 - errmax;
        area = area + area12 - rlist[maxerr];
        if ((defab1 == error1) || (defab2 == error2)) goto _15;
        if ((fabs(rlist[maxerr] - area12) > 1.0e-5 * fabs(area12))
            || (erro12 < .99 * errmax)) goto _10;
        if (extrap) iroff2++;
        else iroff1++;
_10:
        if ((last > 9) && (erro12 > errmax))    /* last > 10 */
            iroff3++;
_15:
        rlist[maxerr] = area1;
        rlist[last] = area2;
        errbnd = max(epsabs,epsrel * fabs(area));

/* Test for roundoff error and eventually set error flag. */
        if (((iroff1 + iroff2) >= 10) || (iroff3 >= 20))
            *ier = 2;
        if (iroff2 > 5)
            ierro = 3;

/* Set error flag in the case that the number of subintervals
    equals limit. */
        if (last == limit)    /* last == limit */
            *ier = 1;

/* Set error flag in the case of bad integrand behavior at some
    points in the integration range. */
        if (max(fabs(a1),fabs(b2)) <= (1.0 +1000.0 * epmach) *
            (fabs(a2) + 1000.0*uflow))
            *ier = 4;

/* Append the newly-created intervals to the list. */
        if (error2 > error1) goto _20;
        alist[last] = a2;
        blist[maxerr] = b1;
        blist[last] = b2;
        elist[maxerr] = error1;
        elist[last] = error2;
        goto _30;
_20:
        alist[maxerr] = a2;
        alist[last] = a1;
        blist[last] = b1;
        rlist[maxerr] = area2;
        rlist[last] = area1;
        elist[maxerr] = error2;
        elist[last] = error1;

/* Call dqsort to maintain the descending ordering in the list of error
    estimates and select the subinterval with nrmax-th largest
    error estimate (to be bisected next). */
_30:
        dqsort(limit,last,&maxerr,&errmax,elist,iord,&nrmax);
        if (errsum <= errbnd) goto _115;
        if (*ier != 0) goto _100;
        if (last == 1) goto _80;    /* last == 2 */
        if (noext) goto _90;        /* goto 90 */
        erlarg -= erlast;
        if (fabs(b1-a1) > small)
            erlarg += erro12;
        if (extrap) goto _40;

/* Test whether the interval to be bisected next is the smallest interval. */
        if ((fabs(blist[maxerr] - alist[maxerr])) > small)
            goto _90;    /* goto 90 */
        extrap = TRUE;
        nrmax = 1;        /* nrmax = 2 */
_40:
        if ((ierro == 3) || (erlarg <= ertest)) goto _60;

/* The smallest interval has the largest error. Before bisecting, decrease
    the sum of the errors over the larger intervals (erlarg) and
        perform extrapolation.) */
        id = nrmax;
        jupbnd = last;
        if (last > (2 + limit/2))
            jupbnd = limit + 3 - last;
        for (k = id;k <= jupbnd; k++) {
            maxerr = iord[nrmax];
            errmax = elist[maxerr];
            if (fabs(blist[maxerr] - alist[maxerr]) > small)
                goto _90;    /* goto 90 */
            nrmax++;
        }

/* Perform extrapolation. */
_60:
        numrl2++;
        rlist2[numrl2] = area;
        reseps=dqext(&numrl2,rlist2,&abseps,res3la,&nres);
        ktmin++;
        if ((ktmin > 5) && (*abserr < 1.0e-3 * errsum)) *ier = 5;
        if (abseps >= *abserr) goto _70;
        ktmin = 0;
        *abserr = abseps;
        result = reseps;
        correc = erlarg;
        ertest = max(epsabs,epsrel * fabs(reseps));
        if (*abserr <= ertest) goto _100;

/* Prepare bisection of the smallest interval. */
_70:
        if (numrl2 == 0) noext = TRUE;
        if (*ier == 5) goto _100;
        maxerr = iord[0];
        errmax = elist[maxerr];
        nrmax = 0;
        extrap = FALSE;
        small = small * 0.5;
        erlarg = errsum;
        goto _90;        /* goto 90 */
_80:
        small = fabs(b-a)*0.375;
        erlarg = errsum;
        ertest = errbnd;
        rlist2[1] = area;
_90:
        ;
    }                    /* 90: */
_100:
    if (*abserr == oflow) goto _115;
    if ((*ier + ierro) == 0) goto _110;
    if (ierro == 3) *abserr += correc;
    if (*ier == 0) *ier = 3;
    if ((result != 0.0) && (area != 0.0)) goto _105;
    if (*abserr > errsum) goto _115;
    if (area == 0.0) goto _130;
    goto _110;
_105:
    if (*abserr/fabs(result) > errsum/fabs(area)) goto _115;

/* Test on divergence. */
_110:
    if ((ksgn == -1) && (max(fabs(result),fabs(area)) <= defabs * .01))
        goto _130;
    if ((0.01 > result/area) || (result/area > 100.0) ||
        (errsum > fabs(area))) *ier = 6;
    goto _130;

/* Compute global integral. */
_115:
    result = 0.0;
    for (k = 0; k <= last; k++)
        result += rlist[k];
    *abserr = errsum;
_130:
    if (*ier > 2) (*ier)--;
_140:
    *neval = 42 * last - 21;
    return result;
}


/*  DQAWC - Computation of Cauchy principal value
 *
 *  PARAMETERS:
 *
 *      f() -   double precision function defining the integrand.
 *
 *      a   -   lower limit of integration
 *
 *      b   -   upper limit of integration
 *
 *      c   -   parameter in the weight function
 *
 *      epsabs  -   absolute accuracy requested
 *
 *      epsrel  -   relative accuracy requested
 *
 *      abserr  -   estimate of the modulus of the absolute error
 *
 *      neval   -   number of function evaluations
 *
 *      ier     -   error code
 */
double dqawc(dq_function_type f,double a,double b,double c,double epsabs,
    double epsrel,double *abserr,int *neval,int *ier, void* user_data)
{
        double result;

        result = dqawce(f,a,b,c,epsabs,epsrel,abserr,neval,ier, user_data);
        return result;
}




/*  DQAWCE - Computation of Cauchy principal value
 *
 *  PARAMETERS:
 *
 *      f() -   double precision function defining the integrand.
 *
 *      a   -   lower limit of integration
 *
 *      b   -   upper limit of integration
 *
 *      c   -   parameter in the weight function
 *
 *      epsabs  -   absolute accuracy requested
 *
 *      epsrel  -   relative accuracy requested
 *
 *      abserr  -   estimate of the modulus of the absolute error
 *
 *      neval   -   number of function evaluations
 *
 *      ier     -   error code
 */
double dqawce(dq_function_type f,double a,double b,double c,double epsabs,
    double epsrel,double *abserr,int *neval,int *ier, void* user_data)
{
    double aa,area,area1,area2,area12,a1,a2,bb,b1,b2;
    double errbnd,errmax,error1,error2,erro12,errsum,result;
    double alist[LIMIT],blist[LIMIT],rlist[LIMIT],elist[LIMIT];

    int iord[LIMIT],iroff1,iroff2,k,krule,last,maxerr,nrmax,nev;
    int limit;

    limit = LIMIT - 1;
    *ier = 6;
    *neval = 0;
    last = 0;
    alist[0] = a;
    blist[0] = b;
    rlist[0] = 0.0;
    elist[0] = 0.0;
    iord[0] = 0;
    result = 0.0;
    *abserr = 0.0;
    if ((c==a) || (c==b) || ((epsabs < 0.0) && (epsrel < 0.0))) goto _999;

/*  First approximation to the integral.    */
    aa = a;
    bb = b;
    if (a <= b) goto _10;
    aa = b;
    bb = a;
_10:
    *ier = 0;
    krule = 1;
    result = dqc25c(f,aa,bb,c,abserr,&krule,neval, user_data);
    rlist[0] = result;
    elist[0] = *abserr;
    iord[0] = 0;
    alist[0] = a;
    blist[0] = b;

/*  Test on accuracy.   */
    errbnd = max(epsabs,epsrel * fabs(result));
    if (limit == 0) *ier = 1;
    if ((*abserr < min(1.0e-2 * fabs(result),errbnd))  || (*ier == 1))
        goto _70;

/*  Initialization. */
    alist[0] = aa;
    blist[0] = bb;
    rlist[0] = result;
    errmax = *abserr;
    maxerr = 0;
    area = result;
    errsum = *abserr;
    nrmax = 0;
    iroff1 = 0;
    iroff2 = 0;

/*  Main loop.  */
    for (last = 1;last < limit;last++) {
/* Bisect the subinterval with nrmax-th largest error estimate.    */
        a1 = alist[maxerr];
        b1 = 0.5 * (alist[maxerr]+blist[maxerr]);
        b2 = blist[maxerr];
        if ((c <= b1) && (c > a1)) b1 = 0.5 * (c + b2);
        if ((c >  b1) && (c < b2)) b1 = 0.5 * (a1 + c);
        a2 = b1;
        krule = 2;
        area1 = dqc25c(f,a1,b1,c,&error1,&krule,&nev, user_data);
        *neval = *neval + nev;
        area2 = dqc25c(f,a2,b2,c,&error2,&krule,&nev, user_data);
        *neval = *neval + nev;

/*  Improve previous approximations to integral and error and
 *  test for accuracy.
 */
        area12 = area1 + area2;
        erro12 = error1 + error2;
        errsum += (erro12-errmax);
        area += (area12-rlist[maxerr]);
        if (fabs(rlist[maxerr]-area12) < (1.0e-5*fabs(area12)) &&
            (erro12 >= 0.99 * errmax) && (krule == 0)) iroff1++;
        if ((last > 10) && (erro12 > errmax) && (krule == 0))
            iroff2++;
        rlist[maxerr] = area1;
        rlist[last] = area2;
        errbnd = max(epsabs,epsrel*fabs(area));
        if (errsum <= errbnd) goto _15;

/*  Test for roundoff error and eventually set error flag. */
        if ((iroff1 >= 6) && (iroff2 > 20)) *ier = 2;

/*  Set error flag in the case that number of interval bisections
 *  exceeds limit.
 */
        if (last == limit) *ier = 1;

/*  Set error flag in the case of bad integrand behavior at a point
 *  of the integration range.
 */
        if (max(fabs(a1),fabs(b2)) <= (1.0 + 1.0e3*epmach)*
            (fabs(a2) + 1.0e3 * uflow)) *ier = 3;
/* Append newly created intervals to list. */
_15:
        if (error2 <= error1) {
            alist[last] = a2;
            blist[maxerr] = b1;
            blist[last] = b2;
            elist[maxerr] = error1;
            elist[last] = error2;
        }
        else {
            alist[maxerr] = a2;
            alist[last] = a1;
            blist[last] = b1;
            rlist[maxerr] = area2;
            rlist[last] = area1;
            elist[maxerr] = error2;
            elist[last] = error1;
        }

/*  Call subroutine dqsort to maintain descending ordering in the list. */
        dqsort(limit,last,&maxerr,&errmax,elist,iord,&nrmax);

/*  Jump out of loop.   */
        if ((*ier != 0) || (errsum <= errbnd)) goto _50;
    }

/*  Compute final result.   */
_50:
    result = 0.0;
    for (k=0;k<=last;k++) {
        result += rlist[k];
    }
    *abserr = errsum;
_70:
    if (aa == b) result = -result;
_999:
    return result;
}




#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>

//#define P   0.9

/* DQAWF - Approximation to Fourier integral. (From QUADPACK)
 *
 *
 * PARAMETERS:
 *
 *    f() - double precision function to be integrated.
 *
 *    a - lower limit of integration.
 *
 *    omega - parameter in weight function.
 *
 *    sincos - indicates which weight function to use:
 *        sincos = COSINE (= 1) --- use cos(omega*x)
 *        sincos = SINE   (= 2) --- use sin(omega*x)
 *
 *    epsabs - absolute accuracy requested.
 *
 */
double dqawf(dq_function_type f,double a,double omega,int sincos,
    double epsabs,double *abserr,int *neval,int *ier, void* user_data)
{
    double **chebmo,erlst[50];
    double result,rslst[50];

    int ierlst[50],i;
    int limlst;

    if ((chebmo = (double **)calloc(MAXP1,sizeof(double *))) == NULL) {
        fprintf(stderr,"Out of memory in dqawf!\n");
        exit(1);
    }
    for (i = 0;i < MAXP1; i++) {
        if ((chebmo[i] = (double *)calloc(25,sizeof(double))) == NULL) {
            fprintf(stderr,"Out of memory in dqawf!\n");
            exit(1);
        }
    }
    *ier = 6;
    *neval = 0;
    result = 0.0;
    *abserr = 0.0;

/* Dimensioning parameters.
 *    limlst - upper bound on number of cycles,
 *    MAXP1 - upper bound on the number of Chebyshev moments.
 */
     limlst = 50;

 /* Check validity of limlst and MAXP1. */
     if ((limlst < 3) || (MAXP1 < 1))
         goto _10;

/* Prepare call for dqawfe. */
    result=dqawfe(f,a,omega,sincos,epsabs,limlst,MAXP1,
        abserr,neval,ier,rslst,erlst,ierlst,chebmo, user_data);
_10:
    for (i = 0; i < MAXP1; i++)
        free(chebmo[i]);
    free(chebmo);
    return result;
}




#define P   0.9

/* DQAWFE - Approximation to Fourier integral. (From QUADPACK)
 *
 *
 * PARAMETERS:
 *
 *    f() - double precision function to be integrated.
 *
 *    a - lower limit of integration.
 *
 *    omega - parameter in weight function.
 *
 *    sincos - indicates which weight function to use:
 *        sincos = COSINE (= 1) --- use cos(omega*x)
 *        sincos = SINE   (= 2) --- use sin(omega*x)
 *
 *    epsabs - absolute accuracy requested.
 *
 *    limlst - upper bound on the number of cycles.
 *
 *  limit - upper bound on the number of subintervals (XXX deleted XXX).
 *
 *    maxp1 - upper bound on the number of Chebyshev moments
 *        which can be stored.
 */
double dqawfe(dq_function_type f,double a,double omega,int sincos,
    double epsabs,int limlst,int maxp1,
    double *abserr,int *neval,int *ier,double rslst[],
    double erlst[],int ierlst[],double **chebmo, void* user_data)
{
    double abseps,correc,cycle,c1,c2,dl,drl;
    double ep,eps,epsa,errsum,fact,p1,psum[52],reseps;
    double res3la[3],result;

    int ktmin,l,ll,momcom,nev,nres,numrl2,lst;

/* Test on validity of parameters. */
    result = 0.0;
    *abserr = 0.0;
    *neval = 0;
    *ier = 0;
    ll = 0;
    if (((sincos != COSINE) && (sincos != SINE)) || (epsabs <= 0.0) ||
        (limlst < 3)) *ier = 6;
    if (*ier == 6) goto _999;
    if (omega != 0.0) goto _10;

/* Integration by DQAGI if omega is zero. */
    if (sincos == COSINE)
        result = dqagi(f,0.0,1,epsabs,0.0,abserr,neval,ier, user_data);
    rslst[0] = result;
    erlst[0] = *abserr;
    ierlst[0] = *ier;
    goto _999;

/* Initialization. */
_10:
    res3la[0] = 0.0;    /* res3la must be initialized to 0.0 */
    res3la[1] = 0.0;
    res3la[2] = 0.0;
    l = fabs(omega);
    dl = 2 * l + 1;
    cycle = dl * Pi / fabs(omega);
    *ier = 0;
    ktmin = 0;
    *neval = 0;
    numrl2 = -1;    /* used as array index. first use is after increment. */
    nres = 0;
    c1 = a;
    c2 = cycle + a;
    p1 = 1.0 - P;
    eps = epsabs;
    if (epsabs > (uflow / p1))
        eps = epsabs * p1;
    ep = eps;
    fact = 1.0;
    correc = 0.0;
    *abserr = 0.0;
    errsum = 0.0;

/* Main Loop */
    for (lst = 0; lst < limlst; lst++) {

/* Integrate over current subinterval. */
/*    dla = lst;  This line is in the original code, but dla is unused. */
    epsa = eps * fact;
    rslst[lst] = dqfour(f,c1,c2,omega,sincos,epsa,0.0,lst+1,maxp1,  // lst+1
        &erlst[lst],&nev,&ierlst[lst],&momcom,chebmo, user_data);
    *neval += nev;
    fact *= P;
    errsum += erlst[lst];
    drl = 50.0 * fabs(rslst[lst]);

/* Test on accuracy with partial sum. */
    if (((errsum + drl) <= epsabs) && (lst >= 5))
        goto _80;
    correc = max(correc,erlst[lst]);
    if (ierlst[lst] != 0)
        eps = max(ep,correc * p1);
    if (ierlst[lst] != 0)
        *ier = 7;
    if ((*ier == 7) && ((errsum + drl) <= (correc * 10.0))
        && (lst > 4)) goto _80;
    numrl2++;
    if (lst > 0)
        goto _20;
    psum[0] = rslst[0];
    goto _40;
_20:
    psum[numrl2] = psum[ll] + rslst[lst];

    if (lst == 1)
        goto _40;

/* Test on maximum number of subintervals. */
    if (lst == limlst-1)
        *ier = 8;

/* Perform new extrapolation. */
    reseps=dqext(&numrl2,psum,&abseps,res3la,&nres);

/* Test whether extrapolated result is influenced by roundoff. */
    ktmin++;
    if ((ktmin >= 15) && (*abserr <= 0.001 * (errsum + drl)))
        *ier = 9;
    if ((abseps > *abserr) && (lst != 2))
        goto _30;
    *abserr = abseps;
    result = reseps;
    ktmin = 0;

/* If ier is not 0, check whether direct result (partial sum) or
 * extrapolated result yields the best integral approximation.
 */
     if (((*abserr + 10.0 * correc) <= epsabs) || (*abserr <= epsabs) &&
         (10.0 * correc >= epsabs)) goto _60;
_30:
    if ((*ier != 0) && (*ier != 7))
        goto _60;
_40:
    ll = numrl2;
    c1 = c2;
    c2 += cycle;
_50:
    ;
    }

 /* Set final result and error estimate. */
_60:
    (*abserr)+=(10.0 * correc);
    if (*ier == 0)
        goto _999;
    if ((result != 0.0) && (psum[numrl2] != 0.0))
        goto _70;
    if (*abserr > errsum)
        goto _80;
    if (psum[numrl2] == 0.0)
        goto _999;
_70:
    if ((*abserr / fabs(result) > (errsum+drl) / fabs(psum[numrl2])))
        goto _80;
    if ((*ier >= 1) && (*ier != 7))
        (*abserr) += drl;
    goto _999;
_80:
    result = psum[numrl2];
    *abserr = errsum + drl;
_999:
    return result;

}


double dqawo(dq_function_type f,double a,double b,double omega, int sincos,
    double epsabs,double epsrel,double *abserr,int *neval,
    int *ier, void* user_data)
{
    double **chebmo,result;
    int i,momcom;

    if ((chebmo = (double **)calloc(MAXP1,sizeof(double *))) == NULL) {
        fprintf(stderr,"Out of memory in DQAWO!\n");
        exit(1);
    }
    for (i = 0;i < MAXP1; i++) {
        if ((chebmo[i] = (double *)calloc(25,sizeof(double))) == NULL) {
            fprintf(stderr,"Out of memory in DQAWO!\n");
            exit(1);
        }
    }

    momcom = 0;
    result = dqfour(f,a,b,omega,sincos,epsabs,epsrel,
        1,MAXP1,abserr,neval,ier,&momcom,chebmo, user_data);
    for (i = 0;i < MAXP1; i++)
        free(chebmo[i]);
    free(chebmo);
    return result;
}


/*  DQAWS - Approximation to integral with algebraic and/or logarithmic
 *          singularities.
 *
 *  PARAMETERS:
 *
 *      f() - double precision function to be integrated.
 *
 *      a   - double lower limit of integration.
 *
 *      b   - upper limit of integration.
 *
 *      alfa - parameter in the weight function.
 *
 *      beta - parameter in the weight function.
 *
 *      wgtfunc - indicates which weight function is to be used.
 *                  = 1:    (x-a)^alfa * (b-x)^beta
 *                  = 2:    (x-a)^alfa * (b-x)^beta * log(x-a)
 *                  = 3:    (x-a)^alfa * (b-x)^beta * log(b-x)
 *                  = 4:    (x-a)^alfa * (b-x)^beta * log(x-a) * log(b-x)
 *
 *      epsabs  - absolute accuracy requested.
 *
 *      epsrel  - relative accuracy requested.
 *
 */
double dqaws(dq_function_type f,double a,double b,double alfa,double beta,
        int wgtfunc,double epsabs,double epsrel,double *abserr,
        int *neval,int *ier, void* user_data)
{
    double result;

    result = dqawse(f,a,b,alfa,beta,wgtfunc,epsabs,epsrel,abserr,
                neval,ier, user_data);
    return result;
}




/*  DQAWSE - Approximation to integral with algebraic and/or logarithmic
 *          singularities.
 *
 *  PARAMETERS:
 *
 *      f() - double precision function to be integrated.
 *
 *      a   - double lower limit of integration.
 *
 *      b   - upper limit of integration.
 *
 *      alfa - parameter in the weight function.
 *
 *      beta - parameter in the weight function.
 *
 *      wgtfunc - indicates which weight function is to be used.
 *                  = 1:    (x-a)^alfa * (b-x)^beta
 *                  = 2:    (x-a)^alfa * (b-x)^beta * log(x-a)
 *                  = 3:    (x-a)^alfa * (b-x)^beta * log(b-x)
 *                  = 4:    (x-a)^alfa * (b-x)^beta * log(x-a) * log(b-x)
 *
 *      epsabs  - absolute accuracy requested.
 *
 *      epsrel  - relative accuracy requested.
 *
 */
double dqawse(dq_function_type f,double a,double b,double alfa,double beta,
        int wgtfunc,double epsabs,double epsrel,double *abserr,
        int *neval,int *ier, void* user_data)
{
    double alist[LIMIT],blist[LIMIT],rlist[LIMIT],elist[LIMIT];
    double ri[25],rj[25],rh[25],rg[25];
    double area,area1,area12,area2,a1,a2,b1,b2,centre;
    double errbnd,errmax,error1,erro12,error2,errsum;
    double resas1,resas2,result;

    int iord[LIMIT],iroff1,iroff2,k,last,limit,maxerr,nev,nrmax;

    limit = LIMIT - 1;
/*  Test on validity of parameters. */
    *ier = 6;
    *neval = 0;
    rlist[0] = 0.0;
    elist[0] = 0.0;
    iord[0] = 0;
    result = 0.0;
    *abserr = 0.0;
    if ((b <= a) || ((epsabs < 0.0) && (epsrel < 0.0)) ||
        (alfa <= -1.0) || (beta <= -1.0) || (wgtfunc < 1) ||
        (wgtfunc > 4) || (limit < 1)) goto _999;
    *ier = 0;

/*  Compute the modified Chebyshev moments. */
    dqmomo(alfa,beta,ri,rj,rg,rh,wgtfunc);

/*  Integrate over the intervals (a,(a+b)/2) and ((a+b)/2,b). */
    centre = 0.5 * (a+b);
    area1 = dqc25s(f,a,b,a,centre,alfa,beta,ri,rj,rg,rh,&error1,
        &resas1,wgtfunc,&nev, user_data);
    *neval = *neval + nev;
    area2 = dqc25s(f,a,b,centre,b,alfa,beta,ri,rj,rg,rh,&error2,
        &resas2,wgtfunc,&nev, user_data);
    *neval = *neval + nev;
    result = area1 + area2;
    *abserr = error1 + error2;

/* Test on accuracy. */
    errbnd = max(epsabs,epsrel * fabs(result));

/*  Initialization. */
    if (error1 >= error2) {
        alist[0] = a;
        alist[1] = centre;
        blist[0] = centre;
        blist[1] = b;
        rlist[0] = area1;
        rlist[1] = area2;
        elist[0] = error1;
        elist[1] = error2;
    }
    else {
        alist[0] = centre;
        alist[1] = a;
        blist[0] = b;
        blist[1] = centre;
        rlist[0] = area2;
        rlist[1] = area1;
        elist[0] = error2;
        elist[1] = error1;
    }
    iord[0] = 0;
    iord[1] = 1;
    if (limit == 1) *ier = 1;
    if ((*abserr <= errbnd) || (*ier == 1)) goto _999;
    errmax = elist[0];
    maxerr = 0;
    nrmax = 0;
    area = result;
    errsum = maxerr;
    iroff1 = 0;
    iroff2 = 0;

/*  Main loop. */
    for (last = 2;last < limit;last++) {
        a1 = alist[maxerr];
        b1 = 0.5 * (alist[maxerr]+blist[maxerr]);
        a2 = b1;
        b2 = blist[maxerr];

        area1 = dqc25s(f,a,b,a1,b1,alfa,beta,ri,rj,rg,rh,&error1,
                &resas1,wgtfunc,&nev, user_data);
        *neval = *neval + nev;
        area2 = dqc25s(f,a,b,a2,b2,alfa,beta,ri,rj,rg,rh,&error2,
                &resas2,wgtfunc,&nev, user_data);
        *neval = *neval + nev;

/*  Improve previous approximation and error test for accuracy. */
        area12 = area1 + area2;
        erro12 = error1 + error2;
        errsum += (erro12 - errmax);
        area += (area12-rlist[maxerr]);
        if ((a == a1) || (b == b2)) goto _30;
        if ((resas1 == error1) || (resas2 == error2)) goto _30;

/*  Test for roundoff error. */
        if ((fabs(rlist[maxerr]-area12) < (1.0e-5 * fabs(area12))) &&
            (erro12 >= (0.99 *errmax))) iroff1++;
        if ((last > 9) && (erro12 > errmax)) iroff2++;
_30:
        rlist[maxerr] = area1;
        rlist[last] = area2;

/*  Test on accuracy. */
        errbnd = max(epsabs,epsrel*fabs(area));
        if (errsum <= errbnd) goto _35;

/*  Set error flag in the case that number of intervals exceeds limit. */
        if (last == limit) *ier = 1;

/*  Set error flag in the case of roundoff error. */
        if ((iroff1 > 5) || (iroff2 > 19)) *ier = 2;

/*  Set error flag in case of bad integrand behavior at interior points. */
        if ( max(fabs(a1),fabs(b2)) <= ((1.0 + 1.0e3 * epmach) *
                (fabs(a2)+1.0e3 * uflow)) ) *ier = 3;
/*  Append the newly created intervals to the list. */
_35:
        if (error2 <= error1) {
            alist[last] = a2;
            blist[maxerr] = b1;
            blist[last] = b2;
            elist[maxerr] = error1;
            elist[last] = error2;
        }
        else {
            alist[maxerr] = a2;
            alist[last] = a1;
            blist[last] = b1;
            rlist[maxerr] = area2;
            rlist[last] = area1;
            elist[maxerr] = error2;
            elist[last] = error1;
        }

/*  Call subroutine qsort to maintain the descending ordering. */
        dqsort(limit,last,&maxerr,&errmax,elist,iord,&nrmax);

/*  Jump out of loop. */
    if ((*ier != 0) || (errsum <= errbnd)) break;
    }
    result = 0.0;
    for (k=0;k<=last;k++) {
        result += rlist[k];
    }
    *abserr = errsum;
_999:
    return result;
}




/*  DQC25C - Integration rules for the computation of Cauchy
 *          principal value integrals.
 *
 *  PARAMETERS:
 *
 *      f() - double precision function defining the integrand
 *
 *      a   - left end point of the integration interval
 *
 *      b   - right end point of the integration interval
 *
 *      c   - parameter in the weight function
 *
 *      abserr  - estimate of the modulus of the absolute error
 *
 *      krul    - key which is decreased by 1 if the 15-point
 *                  Gauss-Kronrod scheme is used
 *
 *      neval   - number of function evaluations
 */
double dqc25c(dq_function_type f,double a,double b,double c,double *abserr,
        int *krul, int *neval, void* user_data)
{
    static double x[11] = {
        0.99144486137381041114,
        0.96592582628906828675,
        0.92387953251128675613,
        0.86602540378443864676,
        0.79335334029123516458,
        0.70710678118654752440,
        0.60876142900872063942,
        0.50000000000000000000,
        0.38268343236508977173,
        0.25881904510252076235,
        0.13052619222005159155};
    double ak22,amom0,amom1,amom2,cc,centr;
    double cheb12[13],cheb24[25],fval[25];
    double hlgth,resabs,resasc,res12,res24,u,result;
    int i,isym,k;
    int unitialized_value = 0xCCCCCCCC;
    int kp = unitialized_value;
    double p2 = unitialized_value, p3 = unitialized_value, p4 = unitialized_value;

    cc = (2.0 * c - b - a) / (b - a);
    if (fabs(cc) < 1.1) goto _10;

/*  Apply the 15-point Gauss-Kronrod scheme.    */
    (*krul)--;
    result = G_K15W(f,dqwgtc,c,p2,p3,p4,kp,a,b,abserr,&resabs,&resasc, user_data);
    *neval = 15;
    if (resasc == *abserr) (*krul)++;
    goto _50;

/*  Use the generalized Clenshaw-Curtis method. */
_10:
    hlgth = 0.5 * (b - a);
    centr = 0.5 * (b + a);
    *neval = 25;
    fval[0] = 0.5 * f(hlgth+centr, user_data);
    fval[12] = f(centr, user_data);
    fval[24] = 0.5 * f(centr-hlgth, user_data);
    for (i=1;i<12;i++) {
        u = hlgth * x[i-1];
        isym = 24 - i;
        fval[i] = f(u+centr, user_data);
        fval[isym] = f(centr-u, user_data);
    }

/*  Compute the Chebyshev series expansion. */
    dqcheb(x,fval,cheb12,cheb24);

/*  The modified Chebyshev moments are computed by forward
 *  recursion, using amom0 and amom1 as starting values.
 */
    amom0 = log(fabs((1.0-cc)/(1.0+cc)));
    amom1 = 2.0 + cc * amom0;
    res12 = cheb12[0] * amom0 + cheb12[1] * amom1;
    res24 = cheb24[0] * amom0 + cheb24[1] * amom1;
    for (k = 2;k < 13;k++) {
        amom2 = 2.0 * cc * amom1 - amom0;
        ak22 = (k-1) * (k-1);
        if ((k/2)*2 != k) amom2 -= (4.0 / (ak22 - 1.0));
        res12 += (cheb12[k] * amom2);
        res24 += (cheb24[k] * amom2);
        amom0 = amom1;
        amom1 = amom2;
    }
    for (k = 13;k < 25;k++) {
        amom2 = 2.0 * cc * amom1 - amom0;
        ak22 = (k-1) * (k-1);
        if ((k/2)*2 != k) amom2 -= (4.0 /(ak22 - 1.0));
        res24 += (cheb24[k] * amom2);
        amom0 = amom1;
        amom1 = amom2;
    }
    result = res24;
    *abserr = fabs(res24-res12);
_50:
    return result;
}




#define NMAC 27

double dqc25o(dq_function_type f,double a,double b,double omega,int sincos,
    int nrmom,int maxp1,int ksave,double *abserr,int *neval,
    double *resabs,double *resasc,int *momcom,double **chebmo, void* user_data)
{
    static double x[11] = {
        0.99144486137381041114,
        0.96592582628906828675,
        0.92387953251128675613,
        0.86602540378443864676,
        0.79335334029123516458,
        0.70710678118654752440,
        0.60876142900872063942,
        0.50000000000000000000,
        0.38268343236508977173,
        0.25881904510252076235,
        0.13052619222005159155};
    double ac,an,an2,as,asap,ass,centr,conc,cons,cospar;
    double estc,ests,hlgth,parint,par2,par22;
    double resc12,resc24,ress12,ress24,result,sinpar;
    double cheb12[13],cheb24[25],d[28],d1[28],d2[28];
    double d3[28],fval[25],v[28];
    int unitialized_value = 0xCCCCCCCC;
    double p2 = unitialized_value, p3 = unitialized_value, p4 = unitialized_value;
    int i,isym,j,k,m,noequ,noeq1,mm1;

    centr = 0.5 * (b + a);
    hlgth = 0.5 * (b - a);
    parint = omega * hlgth;

/* Compute the integral using the 15-point Gauss-Kronrod formula
 * if the value of the parameter in the integrand is small or
 * is less than (bb-aa)/2^(maxp1-2), where (aa,bb) is the original
 * integration interval.
 */
    if (fabs(parint) > 2.0) goto _10;
     result = G_K15W(f,dqwgto,omega,p2,p3,p4,sincos,a,b,
              abserr,resabs,resasc, user_data);
     *neval = 15;
     goto _190;

 /* Compute the integral using the generalized Clenshaw-Curtis method. */
 _10:
    conc = hlgth * cos(centr * omega);
    cons = hlgth * sin(centr * omega);
     *resasc = oflow;
    *neval = 25;

/* Check whether the Chebyshev moments for this interval have
 * already been computed.
 */
     if ((nrmom < *momcom) || (ksave == 1)) goto _140;

/* Compute a new set of Chebyshev moments. */
    m = *momcom + 1;
/*** Add variable mm1 to ease transliteration from FORTRAN array
 *** indexing to C indexing.
 ***/
    mm1 = m - 1;
    par2 = parint * parint;
    par22 = par2 + 2.0;
    sinpar = sin(parint);
    cospar = cos(parint);

/* Compute the Chebyshev moments with respect to cosine. */
    v[0] = 2.0 * sinpar / parint;
    v[1] = (8.0 * cospar + (par2 + par2 - 8.0) * sinpar / parint) / par2;
    v[2] = (32.0 * (par2 - 12.0) * cospar + (2.0 * ((par2 - 80.0) *
        par2 + 192.0) * sinpar) / parint) / (par2 * par2);
    ac = 8.0 * cospar;
    as = 24.0 * parint * sinpar;
    if (fabs(parint) > 24.0) goto _70;

/* Compute the Chebyshev moments as the solutions of a boundary
 * value problem with 1 initial value (v[2]) and 1 end value
 * (computed using an asymptotic formula).
 */

     noequ = NMAC - 3;
     noeq1 = noequ - 1;
     an = 6.0;
    for (k = 0; k <= noeq1; k++) {
         an2 = an * an;
         d[k] = -2.0 * (an2 - 4.0) * (par22 - an2 - an2);
         d2[k] = (an - 1.0) * (an - 2.0) * par2;
         d1[k] = (an + 3.0) * (an + 4.0) * par2;
         v[k+3] = as - (an2 - 4.0) * ac;
         an += 2.0;
     }
     an2 = an * an;
     d[noequ] = -2.0 * (an2 - 4.0) * (par22 - an2 - an2);
     v[noequ+3] = as - (an2 - 4.0) * ac;
     v[3] -= (56.0 * par2 * v[2]);
     ass = parint * sinpar;
    asap = (((((210.0 * par2 -1.0) * cospar - (105.0 * par2 - 63.0) *
         ass) / an2 - (1.0 - 15.0 * par2) * cospar + 15.0 * ass) /
         an2 - cospar + 3.0 * ass) / an2 - cospar) / an2;
     v[noequ+3] -= (2.0 * asap * par2 * (an - 1.0) * (an - 2.0));
/* Solve the tridiagonal system by means of Gaussian elimination
 * with partial pivoting.
 */
     for (i = 0; i <= noequ; i++)
         d3[i] = 0.0;
    d2[noequ] = 0.0;
    for (i = 0; i <= noeq1; i++) {
        if (fabs(d1[i]) <= fabs(d[i])) goto _40;
        an = d1[i];
        d1[i] = d[i];
        d[i] = an;
        an = d2[i];
        d2[i] = d[i+1];
        d[i+1] = an;
        d3[i] = d2[i+1];
        d2[i+1] = 0.0;
        an = v[i+4];
        v[i+4] = v[i+3];
        v[i+3] = an;
_40:
        d[i+1] -= (d2[i] * d1[i] / d[i]);
        d2[i+1] -= (d3[i] * d1[i] / d[i]);
        v[i+4] -= (v[i+3] * d1[i] / d[i]);
    }
    v[noequ+3] /= d[noequ];
    v[noequ+2] = (v[noequ+2] - d2[noeq1] * v[noequ+3]) / d[noeq1];
    for (i = 1; i <= noeq1; i++) {
        k = noequ - i - 1;
        v[k+3] = (v[k+3] - d3[k] * v[k+5] - d2[k] * v[k+4]) / d[k];
    }
    goto _90;

/* Compute the Chebyshev moments by means of forward recursion. */
_70:
    an = 4.0;
    for (i = 3; i < 13; i++) {
        an2 = an * an;
        v[i] = ((an2 - 4.0) * (2.0 * (par22 - an2 - an2) *
            v[i-1] - ac) + as - par2 * (an + 1.0) *
            (an + 2.0) * v[i-2]) / (par2 * (an - 1.0) *
            (an - 2.0));
        an += 2.0;
    }
_90:
    for (j = 0; j < 13; j++)
        chebmo[mm1][2*j] = v[j];

/* Compute the Chebyshev moments with respect to sine. */
    v[0] = 2.0 * (sinpar - parint * cospar) / par2;
    v[1] = (18.0 - 48.0 / par2) * sinpar / par2 +(-2.0 + 48.0 / par2) *
        cospar / parint;
    ac = -24.0 * parint * cospar;
    as = -8.0 * sinpar;
    chebmo[mm1][1] = v[0];
    chebmo[mm1][3] = v[1];
    if (fabs(parint) > 24.0) goto _120;
    for (k = 2; k < 12; k++) {
        an = k+1;
        chebmo[mm1][2*k+1] = - sinpar / (an * (2.0 * an - 2.0)) -
            0.25 * parint * (v[k+1] / an - v[k] / (an - 1.0));
    }
    goto _140;

/* Compute the Chebyshev moments by means of forward recursion. */
_120:
    an = 3.0;
    for (i = 2; i < 12; i++) {
        an2 = an * an;
        v[i] = ((an2 - 4.0) * (2.0 * (par22 - an2 - an2) *
            v[i-1] + as) + ac - par2 * (an + 1.0) *
            (an + 2.0) * v[i-2]) / (par2 * (an - 1.0) *
            (an - 2.0));
        an += 2.0;
        chebmo[mm1][2*i+1] = v[i];
    }
_140:
    if (nrmom < *momcom) {
        m = nrmom + 1;
        mm1 = m - 1;
    }
    if ((*momcom < (maxp1 - 1)) && (nrmom >= (*momcom)))
        (*momcom)++;

/* Compute the coefficients of the Chebyshev expansions of degrees
 * 12 and 24 of the function f.
 */
     fval[0] = 0.5 * f(centr+hlgth, user_data);
     fval[12] = f(centr, user_data);
     fval[24] = 0.5 * f(centr-hlgth, user_data);
     for (i = 1; i < 12; i++) {
         isym = 24 - i;
         fval[i] = f(hlgth*x[i-1]+centr, user_data);
         fval[isym] = f(centr-hlgth*x[i-1], user_data);
     }

     dqcheb(x,fval,cheb12,cheb24);

/* Compute the integral and error estimates. */
    resc12 = cheb12[12] * chebmo[mm1][12];
    ress12 = 0.0;
    estc = fabs(cheb24[24] * chebmo[mm1][24]) + fabs((cheb12[12] -
        cheb24[12]) * chebmo[mm1][12]);
    ests = 0.0;
    k = 10;
    for (j = 0; j < 6; j++) {
        resc12 += (cheb12[k] * chebmo[mm1][k]);
        ress12 += (cheb12[k+1] * chebmo[mm1][k+1]);
        estc += fabs((cheb12[k] - cheb24[k]) * chebmo[mm1][k]);
        ests += fabs((cheb12[k+1] - cheb24[k+1]) * chebmo[mm1][k+1]);
        k -= 2;
    }
    resc24 = cheb24[24] * chebmo[mm1][24];
    ress24 = 0.0;
    *resabs = fabs(cheb24[24]);
    k = 22;
    for (j = 0; j < 12; j++) {
        resc24 += (cheb24[k] * chebmo[mm1][k]);
        ress24 += (cheb24[k+1] * chebmo[mm1][k+1]);
        *resabs += (fabs(cheb24[k]) + fabs(cheb24[k+1]));
        if (j <= 4) {
            estc += (fabs(cheb24[k] * chebmo[mm1][k]));
            ests += (fabs(cheb24[k+1] * chebmo[mm1][k+1]));
        }
        k -= 2;
    }
    *resabs *= fabs(hlgth);
    if (sincos == SINE)
        goto _180;
    result = conc * resc24 - cons * ress24;
    *abserr = fabs(conc * estc) + fabs(cons * ests);
    goto _190;
_180:
    result = conc * ress24 + cons * resc24;
    *abserr = fabs(conc * ests) + fabs(cons * estc);
_190:
     return result;
}




double dqc25s(dq_function_type f,double a,double b,double bl,double br,
    double alfa,double beta,double ri[],double rj[],double rg[],
    double rh[],double *abserr,double *resasc,int wgtfunc,int *nev, void* user_data)
{
    static double x[11] = {
        0.99144486137381041114,
        0.96592582628906828675,
        0.92387953251128675613,
        0.86602540378443864676,
        0.79335334029123516458,
        0.70710678118654752440,
        0.60876142900872063942,
        0.50000000000000000000,
        0.38268343236508977173,
        0.25881904510252076235,
        0.13052619222005159155};
    double centr,dc,factor,fix,hlgth,resabs,res12,res24,u,result;
    double cheb12[13],cheb24[25],fval[25];
    int i,isym;

    *nev = 25;
    if ((bl == a) && ((alfa != 0.0) || (wgtfunc == 2) || (wgtfunc == 4)))
        goto _10;
    if ((br == b) && ((beta != 0.0) || (wgtfunc == 3) || (wgtfunc == 4)))
        goto _140;

/*  If a>bl and b<br, apply the 15-point Gauss-Kronrod scheme. */
    result = G_K15W(f,dqwgts,a,b,alfa,beta,wgtfunc,bl,br,abserr,
                &resabs,resasc, user_data);
    *nev = 15;
    goto _270;

/*  This part is only executed if a = bl.
 *  Compute the Chebyshev series expansion of the following function:
 *  f1 = (0.5*(b+b-br-a)-0.5*(br-a)*x)^beta*
 *          f(0.5*(br-a)*x+0.5*(br+a))
 */
_10:
    hlgth = 0.5 * (br-bl);
    centr = 0.5 * (br+bl);
    fix = b-centr;
    fval[0]  = 0.5 * f(hlgth+centr, user_data)*pow(fix-hlgth,beta);
    fval[12] = f(centr, user_data) * pow(fix,beta);
    fval[24] = 0.5 * f(centr-hlgth, user_data)*pow(fix+hlgth,beta);
    for (i=1;i<12;i++) {
        u = hlgth * x[i-1];
        isym = 24 - i;
        fval[i] = f(u+centr, user_data) * pow(fix-u,beta);
        fval[isym] = f(centr-u, user_data) * pow(fix+u,beta);
    }
    factor = pow(hlgth,alfa+1.0);
    result = 0.0;
    *abserr = 0.0;
    res12 = 0.0;
    res24 = 0.0;
    if (wgtfunc > 2) goto _70;
    dqcheb(x,fval,cheb12,cheb24);

/*  wgtfunc = 1  (or 2) */
    for (i=0;i<13;i++) {
        res12 += (cheb12[i] * ri[i]);
        res24 += (cheb24[i] * ri[i]);
    }
    for (i=13;i<25;i++) {
        res24 += (cheb24[i] * ri[i]);
    }
    if (wgtfunc == 1) goto _130;

/*  wgtfunc = 2 */
    dc = log(br-bl);
    result = res24 * dc;
    res12 = 0.0;
    res24 = 0.0;
    for (i=0;i<13;i++) {
        res12 += (cheb12[i] * rg[i]);
        res24 += (cheb24[i] * rg[i]);
    }
    for (i=13;i<25;i++) {
        res24 += (cheb24[i] * rg[i]);
    }
    goto _130;

/*  Compute the Chebyshev series expansion of the following function:
 *      f4 = f1*log(0.5*(b+b-br-a)-0.5*(br-a)*x)
 */
_70:
    fval[0] *= log(fix-hlgth);
    fval[12] *= log(fix);
    fval[24] *= log(fix+hlgth);
    for (i=1;i<12;i++) {
        u = hlgth * x[i-1];
        isym = 24 - i;
        fval[i] *= log(fix-u);
        fval[isym] *= log(fix+u);
    }
    dqcheb(x,fval,cheb12,cheb24);

/*  wgtfunc = 3  (or 4) */
    for (i=0;i<13;i++) {
        res12 += (cheb12[i] * ri[i]);
        res24 += (cheb24[i] * ri[i]);
    }
    for (i=13;i<25;i++) {
        res24 += (cheb24[i] * ri[i]);
    }
    if (wgtfunc == 3) goto _130;

/*  wgtfunc = 4 */
    dc = log(br-bl);
    result = res24 * dc;
    *abserr = fabs((res24-res12)*dc);
    res12 = 0.0;
    res24 = 0.0;
    for (i=0;i<13;i++) {
        res12 += (cheb12[i] * rg[i]);
        res24 += (cheb24[i] * rg[i]);
    }
    for (i=0;i<13;i++) {
        res24 += (cheb24[i] * rg[i]);
    }
_130:
    result = (result + res24) * factor;
    *abserr = (*abserr + fabs(res24-res12)) * factor;
    goto _270;

/*  This part is executed only if b = br
 *
 *  Compute the Chebyshev series expansion of the following function:
 *
 *  f2 = (0.5 *(b+bl-a-a)+0.5*(b-bl)*x)^alfa *
 *      f(0.5*(b-bl)*x+0.5*(b+bl))
 */
_140:
    hlgth = 0.5 * (b-bl);
    centr = 0.5 * (br+bl);
    fix = centr-a;
    fval[0]  = 0.5 * f(hlgth+centr, user_data) * pow(fix+hlgth,alfa);
    fval[12] = f(centr, user_data) * pow(fix,alfa);
    fval[24] = 0.5 * f(centr-hlgth, user_data) * pow(fix-hlgth,alfa);
    for (i=1;i<12;i++) {
        u = hlgth * x[i-1];
        isym = 24 - i;
        fval[i] = f(u+centr, user_data) * pow(fix+u,alfa);
        fval[isym] = f(centr-u, user_data) * pow(fix-u,alfa);
    }
    factor = pow(hlgth,beta+1.0);
    result = 0.0;
    *abserr = 0.0;
    res12 = 0.0;
    res24 = 0.0;
    if ((wgtfunc == 2) || (wgtfunc == 4)) goto _200;

/*  wgtfunc = 1  (or 3)  */
    dqcheb(x,fval,cheb12,cheb24);
    for (i=0;i<13;i++) {
        res12 += (cheb12[i] * rj[i]);
        res24 += (cheb24[i] * rj[i]);
    }
    for (i=13;i<25;i++) {
        res24 += (cheb24[i] * rj[i]);
    }
    if (wgtfunc == 1) goto _260;

/*  wgtfunc = 3  */
    dc = log(br-bl);
    result = res24 * dc;
    *abserr = fabs((res24 - res12) * dc);
    res12 = 0.0;
    res24 = 0.0;
    for (i=0;i<13;i++) {
        res12 += (cheb12[i] * rh[i]);
        res24 += (cheb24[i] * rh[i]);
    }
    for (i=13;i<25;i++) {
        res24 += (cheb24[i] * rh[i]);
    }
_190:
    goto _260;

/*  Compute the Chebyshev series expansion of the following function:
 *
 *      f3 = f2 * log(0.5*(b-bl)*x+0.5*(b+bl-a-a))
 */
_200:
    fval[0] *= log(hlgth+fix);
    fval[12] *= log(fix);
    fval[24] *= log(fix-hlgth);
    for (i=1;i<12;i++) {
        u = hlgth * x[i-1];
        isym = 24 - i;
        fval[i] *= log(u+fix);
        fval[isym] *= log(fix-u);
    }
    dqcheb(x,fval,cheb12,cheb24);

/*  wgtfunc = 2  (or 4)  */
    for (i=0;i<13;i++) {
        res12 += (cheb12[i] * rj[i]);
        res24 += (cheb24[i] * rj[i]);
    }
    for (i=13;i<25;i++) {
        res24 += (cheb24[i] * rj[i]);
    }
    if (wgtfunc == 2) goto _260;
    dc = log(br-bl);
    result = res24 * dc;
    *abserr = fabs((res24-res12) * dc);
    res12 = 0.0;
    res24 = 0.0;

/*  wgtfunc == 4  */
    for (i=0;i<13;i++) {
        res12 += (cheb12[i] * rh[i]);
        res24 += (cheb24[i] * rh[i]);
    }
    for (i=13;i<25;i++) {
        res24 += (cheb24[i] * rh[i]);
    }
_260:
    result = (result + res24)* factor;
    *abserr = (*abserr + fabs(res24-res12))*factor;
_270:
     return result;
}


//#include "cquadpack.h"

void dqcheb(double x[],double fval[],double cheb12[],double cheb24[])
{
    double alam,alam1,alam2,part1,part2,part3;
    double v[12];
    int i,j;

/* Dimensions of input vectors are:
 *        x[11], fval[25], cheb12[13], cheb24[25]
 */
     for (i = 0; i < 12; i++) {
         j = 24 - i;
         v[i] = fval[i] - fval[j];
         fval[i] += fval[j];
     }
     alam1 = v[0] - v[8];
     alam2 = x[5] * (v[2] - v[6] - v[10]);
     cheb12[3] = alam1 + alam2;
     cheb12[9] = alam1 - alam2;
     alam1 = v[1] - v[7] - v[9];
     alam2 = v[3] - v[5] - v[11];
     alam = x[2] * alam1 + x[8] * alam2;
     cheb24[3] = cheb12[3] + alam;
     cheb24[21] = cheb12[3] - alam;
     alam = x[8] * alam1 - x[2] * alam2;
     cheb24[9] = cheb12[9] + alam;
     cheb24[15] = cheb12[9] - alam;
     part1 = x[3] * v[4];
     part2 = x[7] * v[8];
     part3 = x[5] * v[6];
     alam1 = v[0] + part1 + part2;
     alam2 = x[1] * v[2] + part3 + x[9] * v[10];
     cheb12[1] = alam1 + alam2;
     cheb12[11] = alam1 - alam2;
     alam = x[0] * v[1] + x[2] * v[3] + x[4] * v[5] +
         x[6] * v[7] + x[8] * v[9] + x[10] * v[11];
     cheb24[1] = cheb12[1] + alam;
     cheb24[23] = cheb12[1] - alam;
     alam = x[10] * v[1] - x[8] * v[3] + x[6] * v[5] -
         x[4] * v[7] + x[2] * v[9] - x[0] * v[11];
     cheb24[11] = cheb12[11] + alam;
     cheb24[13] = cheb12[11] - alam;
     alam1 = v[0] - part1 + part2;
     alam2 = x[9] * v[2] - part3 + x[1] * v[10];
     cheb12[5] = alam1 + alam2;
     cheb12[7] = alam1 - alam2;
     alam = x[4] * v[1] - x[8] * v[3] - x[0] * v[5] -
         x[10] * v[7] + x[2] * v[9] + x[6] * v[11];
     cheb24[5] = cheb12[5] + alam;
     cheb24[19] = cheb12[5] - alam;
     alam = x[6] * v[1] - x[2] * v[3] - x[10] * v[5] +
         x[0] * v[7] - x[8] * v[9] - x[4] * v[11];
     cheb24[7] = cheb12[7] + alam;
     cheb24[17] = cheb12[7] - alam;
     for (i = 0; i < 6; i++) {
         j = 12 - i;
         v[i] = fval[i] - fval[j];
         fval[i] += fval[j];
     }
     alam1 = v[0] + x[7] * v[4];
     alam2 = x[3] * v[2];
     cheb12[2] = alam1 + alam2;
     cheb12[10] = alam1 - alam2;
     cheb12[6] = v[0] - v[4];
     alam = x[1] * v[1] + x[5] * v[3] + x[9] * v[5];
     cheb24[2] = cheb12[2] + alam;
     cheb24[22] = cheb12[2] - alam;
     alam = x[5] * (v[1] - v[3] - v[5]);
     cheb24[6] = cheb12[6] + alam;
     cheb24[18] = cheb12[6] - alam;
     alam = x[9] * v[1] - x[5] * v[3] + x[1] * v[5];
     cheb24[10] = cheb12[10] + alam;
     cheb24[14] = cheb12[10] - alam;
     for (i = 0; i < 3; i++) {
         j = 6 - i;
         v[i] = fval[i] -fval[j];
         fval[i] += fval[j];
     }
     cheb12[4] = v[0] + x[7] * v[2];
     cheb12[8] = fval[0] - x[7] * fval[2];
     alam = x[3] * v[1];
     cheb24[4] = cheb12[4] + alam;
     cheb24[20] = cheb12[4] - alam;
     alam = x[7] * fval[1] - fval[3];
     cheb24[8] = cheb12[8] + alam;
     cheb24[16] = cheb12[8] - alam;
     cheb12[0] = fval[0] + fval[2];
     alam = fval[1] + fval[3];
     cheb24[0] = cheb12[0] + alam;
     cheb24[24] = cheb12[0] - alam;
     cheb12[12] = v[0] - v[2];
     cheb24[12] = cheb12[12];
      alam = 1.0 / 6.0;
      for (i = 1; i < 12; i++)
          cheb12[i] *= alam;
      alam *= 0.5;
      cheb12[0] *= alam;
      cheb12[12] *= alam;
      for (i = 1; i < 24; i ++)
          cheb24[i] *= alam;
      cheb24[0] *= (0.5 * alam);
      cheb24[24] *= (0.5 * alam);
}





double dqext(int *n,double epstab[],double *abserr,
    double res3la[],int *nres)
{
    double delta1,delta2,delta3,epsinf;
    double error,err1,err2,err3,e0,e1,e1abs,e2,e3;
    double res,result,ss,tol1,tol2,tol3;
    int NN,i,ib,ib2,ie,indx,k1,k2,k3,limexp,newelm,num;

    (*nres)++;
    NN = *n;
    NN++;   /* make NN a FORTRAN array index */
    *abserr = oflow;
    result = epstab[*n];
    if (NN < 3) goto _100;        /* N < 3 */
    limexp = 50;            /* limexp = 50 */
    epstab[*n+2] = epstab[*n];
    newelm = (*n)/2;      /* (n-1)/2 */
    epstab[*n] = oflow;
    num = NN;
    k1 = NN;
    for (i = 1; i <= newelm; i++) {
        k2 = k1 - 1;
        k3 = k1 - 2;
        res = epstab[k1+1];
        e0 = epstab[k3-1];
        e1 = epstab[k2-1];
        e2 = res;
        e1abs = fabs(e1);
        delta2 = e2 - e1;
        err2 = fabs(delta2);
        tol2 = max(fabs(e2),e1abs) * epmach;
        delta3 = e1 - e0;
        err3 = fabs(delta3);
        tol3 = max(e1abs,fabs(e0)) * epmach;
        if ((err2 > tol2) || (err3 > tol3)) goto _10;
        result = res;
        *abserr = err2 + err3;
        goto _100;
_10:
        e3 = epstab[k1-1];
        epstab[k1-1] = e1;
        delta1 = e1 - e3;
        err1 = fabs(delta1);
        tol1 = max(e1abs,fabs(e3)) * epmach;
        if ((err1 <= tol1) || (err2 <= tol2) || (err3 <= tol3)) goto _20;
        ss = 1.0/delta1 + 1.0/delta2 - 1.0/delta3;
        epsinf = fabs(ss*e1);
        if (epsinf > 1.0e-4) goto _30;
_20:
        NN = i + i - 1;
        goto _50;
_30:
        res = e1 + 1.0 / ss;
        epstab[k1-1] = res;
        k1 -= 2;
        error = err2 + fabs(res - e2) + err3;
        if (error > (*abserr)) goto _40;
        *abserr = error;
        result = res;
_40:
        ;
    }
_50:
    if (NN == limexp) NN = 2 * (limexp/2) - 1;
    ib = 1;                        /* ib = 1 */
    if (((num/2) * 2 ) == num) ib = 2;        /* ib = 2 */
    ie = newelm + 1;
    for (i = 1;i <= ie; i++) {
        ib2 = ib + 2;
        epstab[ib-1] = epstab[ib2-1];
        ib = ib2;
    }
    if (num == NN) goto _80;
    indx = num - NN + 1;
    for (i = 1;i <= NN; i++) {
        epstab[i-1] = epstab[indx-1];
        indx++;
    }
_80:
    if (*nres > 3) goto _90;       /* nres >= 4 */
    res3la[(*nres)-1] = result;
    *abserr = oflow;
    goto _100;
_90:
    *abserr = fabs(result - res3la[2]) + fabs(result - res3la[1]) +
        fabs(result - res3la[0]);
    res3la[0] = res3la[1];
    res3la[1] = res3la[2];
    res3la[2] = result;
_100:
    *abserr = max(*abserr,5.0 * epmach * fabs(result));
    *n = NN - 1;
    return result;
    }



/* DQFOUR - Computation of oscillatory integrals.
 *
 *    Calculates an approximation to a given definite integral.
 *        I = integral of F(X) * W(X) over (A,B)
 *            where W(X) = COS(OMEGA * X)
 *               or W(X) = SIN(OMEGA * X)
 *
 * PARAMETERS:
 *
 *    f() - double precision function to be integrated.
 *
 *    a - lower limit of integration.
 *
 *    b - upper limit of integration.
 *
 *    omega - parameter in the integrand weight function.
 *
 *    sincos - indicates which weight function to use:
 *        sincos = COSINE,    W(X) = COS(OMEGA * X)
 *        sincos = SINE,        W(X) = SIN(OMEGA * X)
 *
 *    epsabs - absolute accuracy requested.
 *
 *    epsrel - relative accuracy requested.
 *
 *    limit - upper bound on number of of subdivisions in
 *        the partition of (A,B)
 *
 */

double dqfour(dq_function_type f,double a,double b,double omega,
    int sincos,double epsabs,double epsrel,
    int icall,int maxp1,double *abserr,
    int *neval,int *ier,
    int *momcom,double **chebmo, void* user_data)
{
    double abseps,area,area1,area12,area2;
    double a1,a2,b1,b2,correc,defabs,defab1;
    double defab2,domega,dres,erlarg,erlast,errbnd;
    double errmax,error1,error2,erro12,errsum,ertest;
    double resabs,reseps,result,res3la[3];
    double alist[LIMIT],blist[LIMIT],rlist[LIMIT];
    double elist[LIMIT],rlist2[52],small,width;

    int id,ierro,iroff1,iroff2,iroff3,jupbnd,k,ksgn,limit;
    int ktmin,last,maxerr,nev,nres,nrmax,nrmom,numrl2;
    int extrap,noext,extall,iord[LIMIT],nnlog[LIMIT];

    limit = LIMIT - 1;
/* Test validity of parameters. */
    *ier = 0;
    *neval = 0;
//    last = 0;
    result = 0.0;
    *abserr = 0.0;
    alist[0] = a;
    blist[0] = b;
    rlist[0] = 0.0;
    elist[0] = 0.0;
    iord[0] = 0;
    nnlog[0] = 0;
    if (((sincos != COSINE) && (sincos != SINE)) || ((epsabs < 0.0) &&
        (epsrel < 0.0)) || (icall < 1) || (maxp1 < 1)) *ier = 6;
    if (*ier == 6) return result;

/* First approximation to the integral. */
    domega = fabs(omega);
    nrmom = 0;
    if (icall <= 1)
        *momcom = 0;
_5:
    result = dqc25o(f,a,b,domega,sincos,nrmom,maxp1,0,
        abserr,neval,&defabs,&resabs,momcom,chebmo, user_data);
/* Test on accuracy. */
    dres = fabs(result);
    errbnd = max(epsabs,epsrel*dres);
    rlist[0] = result;
    elist[0] = *abserr;
    iord[0] = 0;
    if ((*abserr <= 100.0 * epmach * defabs) && (*abserr > errbnd))
        *ier = 2;
    if (limit == 0) *ier = 1;
    if ((*ier != 0) || (*abserr <= errbnd) && (*abserr != resabs) ||
        (*abserr == 0.0)) goto _200;

/* Initialization. */
    errmax = *abserr;
    maxerr = 0;             /* maxerr = 1 */
    area = result;
    errsum = *abserr;
    *abserr = oflow;
    nrmax = 0;
    extrap = FALSE;
    noext = FALSE;
    ierro = 0;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ktmin = 0;
    small = fabs(b-a) * 0.75;
    nres = 0;
    numrl2 = -1;
    extall = FALSE;
    if ((0.5 * fabs(b-a) * domega) > 2.0)
        goto _10;
    numrl2 = 0;
    extall = TRUE;
    rlist2[0] = result;
_10:
    if ((0.25 * fabs(b-a) * domega) <= 2.0)
        extall = TRUE;
    ksgn = -1;
    if (dres > (1.0 - 50.0 * epmach) * defabs)
        ksgn = 1;

/* Main loop. */
    for (last = 1; last < limit; last++) {

/* Bisect the interval with the nrmax-th largest error estimate. */
        nrmom = nnlog[maxerr] + 1;
        a1 = alist[maxerr];
        b1 = 0.5 * (alist[maxerr] + blist[maxerr]);
        a2 = b1;
        b2 = blist[maxerr];
        erlast = errmax;
        area1 = dqc25o(f,a1,b1,domega,sincos,nrmom,maxp1,0,
            &error1,&nev,&resabs,&defab1,momcom,chebmo, user_data);
        *neval += nev;
        area2 = dqc25o(f,a2,b2,domega,sincos,nrmom,maxp1,1,
            &error2,&nev,&resabs,&defab2,momcom,chebmo, user_data);
        *neval += nev;

/* Improve previous approximations to integral and error
      and test for accuracy. */
        area12 = area1 + area2;
        erro12 = error1 + error2;
        errsum = errsum + erro12 - errmax;
        area = area + area12 - rlist[maxerr];
        if ((defab1 == error1) || (defab2 == error2)) goto _25;
        if ((fabs(rlist[maxerr] - area12) > 1.0e-5 * fabs(area12))
            || (erro12 < .99 * errmax)) goto _20;
        if (extrap) iroff2++;
        else iroff1++;
_20:
        if ((last > 9) && (erro12 > errmax))    /* last > 10 */
            iroff3++;
_25:
        rlist[maxerr] = area1;
        rlist[last] = area2;
        nnlog[maxerr] = nrmom;
        nnlog[last] = nrmom;
        errbnd = max(epsabs,epsrel * fabs(area));

/* Test for roundoff error and eventually set error flag. */
        if (((iroff1 + iroff2) >= 10) || (iroff3 >= 20))
            *ier = 2;
        if (iroff2 >= 5)
            *ier = 3;

/* Set error flag in the case that the number of subintervals
    equals limit. */
        if (last == limit)    /* last == limit */
            *ier = 1;

/* Set error flag in the case of bad integrand behavior at some
    points in the integration range. */
        if (max(fabs(a1),fabs(b2)) <= (1.0 +1000.0 * epmach) *
            (fabs(a2) + 1000.0*uflow))
            *ier = 4;

/* Append the newly-created intervals to the list. */
        if (error2 <= error1) {
            alist[last] = a2;
            blist[maxerr] = b1;
            blist[last] = b2;
            elist[maxerr] = error1;
            elist[last] = error2;
        }
        else {
            alist[maxerr] = a2;
            alist[last] = a1;
            blist[last] = b1;
            rlist[maxerr] = area2;
            rlist[last] = area1;
            elist[maxerr] = error2;
            elist[last] = error1;
        }
/* Call dqsort to maintain the descending ordering in the list of error
    estimates and select the subinterval with nrmax-th largest
    error estimate (to be bisected next). */

        dqsort(limit,last,&maxerr,&errmax,elist,iord,&nrmax);
        if (errsum <= errbnd) goto _170;
        if (*ier != 0) goto _150;
        if ((last == 1) && (extall)) goto _120;    /* last == 2 */
        if (noext) goto _140;
        if (!extall) goto _50;
        erlarg -= erlast;
        if (fabs(b1-a1) > small)
            erlarg += erro12;
        if (extrap) goto _70;

/* Test whether the interval to be bisected next is the smallest interval. */
_50:
        width = fabs(blist[maxerr] - alist[maxerr]);
        if (width > small)
            goto _140;
        if (extall)
            goto _60;

/* Test whether we can start with the extrapolation procedure (we do
 * this if we integrate over the next interval with use of a Gauss-
 * Kronrod rule) - see routine dqc25o. */
         small *= 0.5;
         if ((0.25 * width * domega) > 2.0)
             goto _140;
         extall = TRUE;
         goto _130;
_60:
        extrap = TRUE;
        nrmax = 1;        /* FORTRAN: nrmax = 2 */
_70:
        if ((ierro == 3) || (erlarg <= ertest))
            goto _90;

/* The smallest interval has the largest error. Before bisecting, decrease
    the sum of the erorrs over the larger intervals (erlarg) and
        perform extrapolation. */
        jupbnd = last;
        if (last > (2 + limit/2))
            jupbnd = limit + 3 - last;
        id = nrmax;
        for (k = id;k <= jupbnd; k++) {
            maxerr = iord[nrmax];
            errmax = elist[maxerr];
            if (fabs(blist[maxerr] - alist[maxerr]) > small)
                goto _140;
            nrmax++;
        }

/* Perform extrapolation. */
_90:
        numrl2++;
        rlist2[numrl2] = area;
        if (numrl2 < 2)
            goto _110;
        reseps = dqext(&numrl2,rlist2,&abseps,res3la,&nres);
        ktmin++;
        if ((ktmin > 5) && (*abserr < 1.0e-3 * errsum)) *ier = 5;
        if (abseps >= *abserr) goto _100;
        ktmin = 0;
        *abserr = abseps;
        result = reseps;
        correc = erlarg;
        ertest = max(epsabs,epsrel * fabs(reseps));
        if (*abserr <= ertest) goto _150;

/* Prepare bisection of the smallest interval. */
_100:
        if (numrl2 == 0) noext = TRUE;
        if (*ier == 5) goto _150;
_110:
        maxerr = iord[0];
        errmax = elist[maxerr];
        nrmax = 0;
        extrap = FALSE;
        small *= 0.5;
        erlarg = errsum;
        goto _140;
_120:
        small *= 0.5;
        numrl2++;
        rlist2[numrl2] = area;
_130:
        erlarg = errsum;
        ertest = errbnd;
_140:
        ;
    }

/* Set the final result. */
_150:
    if ((*abserr == oflow) || (nres == 0)) goto _170;
    if ((*ier + ierro) == 0) goto _165;
    if (ierro == 3) *abserr += correc;
    if (*ier == 0) *ier = 3;
    if ((result != 0.0) && (area != 0.0)) goto _160;
    if (*abserr > errsum) goto _170;
    if (area == 0.0) goto _190;
    goto _165;
_160:
    if (*abserr/fabs(result) > errsum/fabs(area)) goto _170;

/* Test on divergence. */
_165:
    if ((ksgn == -1) && (max(fabs(result),fabs(area)) <= defabs * .01))
        goto _190;
    if ((0.01 > result/area) || (result/area > 100.0) ||
        (errsum > fabs(area))) *ier = 6;
    goto _190;

/* Compute global integral. */
_170:
    result = 0.0;
    for (k = 0; k <= last; k++)
        result += rlist[k];
    *abserr = errsum;
_190:
    if (*ier > 2) (*ier)--;
_200:
    if ((sincos == SINE) && (omega < 0.0))
        result = - result;
    return result;
}




double G_K15(dq_function_type f,double a,double b,double *abserr,
    double *resabs,double *resasc, void* user_data)
{
    static double XGK15[8] = {
        0.99145537112081263921,
        0.94910791234275852453,
        0.86486442335976907279,
        0.74153118559939443986,
        0.58608723546769113029,
        0.40584515137739716691,
        0.20778495500789846760,
        0.00000000000000000000};
    static double WGK15[8] = {
        0.02293532201052922496,
        0.06309209262997855329,
        0.10479001032225018384,
        0.14065325971552591875,
        0.16900472663926790283,
        0.19035057806478540991,
        0.20443294007529889241,
        0.20948214108472782801};
    static double WG7[4] = {
        0.12948496616886969327,
        0.27970539148927666790,
        0.38183005050511894495,
        0.41795918367346938776};
    double fv1[7],fv2[7];
    double absc,centr,dhlgth;
    double fc,fsum,fval1,fval2,hlgth;
    double resg,resk,reskh,result;
    int j,jtw,jtwm1;

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    fc=(*f)(centr, user_data);
    resg = fc * WG7[3];
    resk = fc * WGK15[7];
    *resabs = fabs(resk);
    for (j = 0; j < 3; j++) {
        jtw = 2 * j + 1;
        absc = hlgth * XGK15[jtw];
        fval1 = (*f)(centr-absc, user_data);
        fval2 = (*f)(centr + absc, user_data);
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg += WG7[j] * fsum;
        resk += WGK15[jtw] * fsum;
        *resabs = *resabs + WGK15[jtw] * (fabs(fval1) + fabs(fval2));
    }
    for (j = 0; j < 4; j++) {
        jtwm1 = j * 2;
        absc = hlgth * XGK15[jtwm1];
        fval1 = (*f)(centr-absc, user_data);
        fval2 = (*f)(centr+absc, user_data);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk = resk + WGK15[jtwm1] * fsum;
        *resabs = (*resabs) + WGK15[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    reskh = resk * 0.5;
    *resasc = WGK15[7] * fabs(fc - reskh);
    for (j = 0; j < 7; j++ )
        *resasc = (*resasc) + WGK15[j] * (fabs(fv1[j] - reskh) +
            fabs(fv2[j] - reskh));
    result = resk * hlgth;
    *resabs = (*resabs) * dhlgth;
    *resasc = (*resasc) * dhlgth;
    *abserr = fabs((resk - resg) * hlgth);
    if ((*resasc != 0.0) && (*abserr != 0.0))
        *abserr = (*resasc) * min(1.0,pow((200.0 * (*abserr)/(*resasc)),1.5));
    if (*resabs > uflow/(50.0 * epmach))
        *abserr = max(epmach * 50.0 * (*resabs),(*abserr));
    return result;
}




double G_K15I(dq_function_type f, double boun, int inf, double a, double b,
    double *abserr,double *resabs,double *resasc, void* user_data)
{
    static double XGK15[8] = {
        0.99145537112081263921,
        0.94910791234275852453,
        0.86486442335976907279,
        0.74153118559939443986,
        0.58608723546769113029,
        0.40584515137739716691,
        0.20778495500789846760,
        0.00000000000000000000};
    static double WGK15[8] = {
        0.02293532201052922496,
        0.06309209262997855329,
        0.10479001032225018384,
        0.14065325971552591875,
        0.16900472663926790283,
        0.19035057806478540991,
        0.20443294007529889241,
        0.20948214108472782801};
    static double WG7[4] = {
        0.12948496616886969327,
        0.27970539148927666790,
        0.38183005050511894495,
        0.41795918367346938776};
    double fv1[8],fv2[8];
    double absc,absc1,absc2,centr,dinf;
    double fc,fsum,fval1,fval2,hlgth,resg,resk;
    double reskh,result,tabsc1,tabsc2;
    int j;

    dinf = min((double)(1.0),(double)inf);
    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    tabsc1 = boun + dinf * (1.0 - centr)/centr;
    fval1 = (*f)(tabsc1, user_data);
    if (inf == 2)
        fval1 += (*f)(-tabsc1, user_data);
    fc=(fval1/centr)/centr;
    resg = fc * WG7[3];
    resk = fc * WGK15[7];
    *resabs = fabs(resk);
    for (j = 0; j < 7; j++) {
        absc = hlgth * XGK15[j];
        absc1 = centr - absc;
        absc2 = centr + absc;
        tabsc1 = boun + dinf * (1.0 - absc1)/absc1;
        tabsc2 = boun + dinf * (1.0 - absc2)/absc2;
        fval1 = (*f)(tabsc1, user_data);
        fval2 = (*f)(tabsc2, user_data);
        if (inf == 2) {
            fval1 += (*f)(-tabsc1, user_data);
            fval2 += (*f)(-tabsc2, user_data);
        }
        fval1 = (fval1/absc1)/absc1;
        fval2 = (fval2/absc2)/absc2;
        fv1[j] = fval1;
        fv2[j] = fval2;
        fsum = fval1 + fval2;
        if (j & 1) resg += WG7[j/2] * fsum; /* odd 'j's are truncated */
        resk += WGK15[j] * fsum;
        *resabs = (*resabs) + WGK15[j] * (fabs(fval1) + fabs(fval2));
    }
    reskh = resk * 0.5;
    *resasc = WGK15[7] * fabs(fc - reskh);
    for (j = 0; j < 7; j++ )
        *resasc = (*resasc) + WGK15[j] * (fabs(fv1[j] - reskh) +
            fabs(fv2[j] - reskh));
    result = resk * hlgth;
    *resabs = (*resabs) * hlgth;
    *resasc = (*resasc) * hlgth;
    *abserr = fabs((resk - resg) * hlgth);
    if ((*resasc != 0.0) && (*abserr != 0.0))
        *abserr = (*resasc) * min(1.0,pow((200.0 * (*abserr)/(*resasc)),1.5));
    if (*resabs > uflow/(50.0 * epmach))
        *abserr = max(epmach * 50.0 * (*resabs),(*abserr));
    return result;
}


/*
   originally and errorenously the second argument was declared folowing way: double(*w)()
*/

double G_K15W(dq_function_type f,double (*w)(double, double, double, double, double, int),double p1,double p2,double p3,
    double p4,int kp,double a,double b,double *abserr,
    double *resabs,double *resasc, void* user_data)
{
    static double XGK15[8] = {
        0.99145537112081263921,
        0.94910791234275852453,
        0.86486442335976907279,
        0.74153118559939443986,
        0.58608723546769113029,
        0.40584515137739716691,
        0.20778495500789846760,
        0.00000000000000000000};
    static double WGK15[8] = {
        0.02293532201052922496,
        0.06309209262997855329,
        0.10479001032225018384,
        0.14065325971552591875,
        0.16900472663926790283,
        0.19035057806478540991,
        0.20443294007529889241,
        0.20948214108472782801};
    static double WG7[4] = {
        0.12948496616886969327,
        0.27970539148927666790,
        0.38183005050511894495,
        0.41795918367346938776};
    double fv1[7],fv2[7];
    double absc,absc1,absc2,centr,dhlgth;
    double fc,fsum,fval1,fval2,hlgth;
    double resg,resk,reskh,result;
    int j,jtw,jtwm1;

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    fc=(*f)(centr, user_data) * (*w)(centr,p1,p2,p3,p4,kp);
    resg = fc * WG7[3];
    resk = fc * WGK15[7];
    *resabs = fabs(resk);
    for (j = 0; j < 3; j++) {
        jtw = 2 * j + 1;
        absc = hlgth * XGK15[jtw];
        absc1 = centr - absc;
        absc2 = centr + absc;
        fval1 = (*f)(absc1, user_data) * (*w)(absc1, p1, p2, p3, p4, kp);
        fval2 = (*f)(absc2, user_data) * (*w)(absc2, p1, p2, p3, p4, kp);
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg += WG7[j] * fsum;
        resk += WGK15[jtw] * fsum;
        *resabs = *resabs + WGK15[jtw] * (fabs(fval1) + fabs(fval2));
    }
    for (j = 0; j < 4; j++) {
        jtwm1 = j * 2;
        absc = hlgth * XGK15[jtwm1];
        absc1 = centr - absc;
        absc2 = centr + absc;
        fval1 = (*f)(absc1, user_data) * (*w)(absc1, p1, p2, p3, p4, kp);
        fval2 = (*f)(absc2, user_data) * (*w)(absc2, p1, p2, p3, p4, kp);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk = resk + WGK15[jtwm1] * fsum;
        *resabs = (*resabs) + WGK15[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    reskh = resk * 0.5;
    *resasc = WGK15[7] * fabs(fc - reskh);
    for (j = 0; j < 7; j++ )
        *resasc = (*resasc) + WGK15[j] * (fabs(fv1[j] - reskh) +
            fabs(fv2[j] - reskh));
    result = resk * hlgth;
    *resabs = (*resabs) * dhlgth;
    *resasc = (*resasc) * dhlgth;
    *abserr = fabs((resk - resg) * hlgth);
    if ((*resasc != 0.0) && (*abserr != 0.0))
        *abserr = (*resasc) * min(1.0,pow((200.0 * (*abserr)/(*resasc)),1.5));
    if (*resabs > uflow/(50.0 * epmach))
        *abserr = max(epmach * 50.0 * (*resabs),(*abserr));
    return result;
}




double G_K21(dq_function_type f,double a,double b,double *abserr,
    double *resabs,double *resasc, void* user_data)
{
    static double XGK21[11] = {
        0.99565716302580808074,
        0.97390652851717172008,
        0.93015749135570822600,
        0.86506336668898451073,
        0.78081772658641689706,
        0.67940956829902440623,
        0.56275713466860468334,
        0.43339539412924719080,
        0.29439286270146019813,
        0.14887433898163121088,
        0.00000000000000000000};
    static double WGK21[11] = {
        0.01169463886737187428,
        0.03255816230796472748,
        0.05475589657435199603,
        0.07503967481091995277,
        0.09312545458369760554,
        0.10938715880229764190,
        0.12349197626206585108,
        0.13470921731147332593,
        0.14277593857706008080,
        0.14773910490133849137,
        0.14944555400291690566};
    static double WG10[5] = {
        0.06667134430868813759,
        0.14945134915058059315,
        0.21908636251598204400,
        0.26926671930999635509,
        0.29552422471475287017};
    double fv1[10],fv2[10];
    double absc,centr,dhlgth;
    double fc,fsum,fval1,fval2,hlgth;
    double resg,resk,reskh,result;
    int j,jtw,jtwm1;

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    resg = 0.0;
    fc=(*f)(centr, user_data);
    resk = fc * WGK21[10];
    *resabs = fabs(resk);
    for (j = 0; j < 5; j++) {
        jtw = 2 * j + 1;
        absc = hlgth * XGK21[jtw];
        fval1 = (*f)(centr-absc, user_data);
        fval2 = (*f)(centr+absc, user_data);
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg += WG10[j] * fsum;
        resk += WGK21[jtw] * fsum;
        *resabs = *resabs + WGK21[jtw] * (fabs(fval1) + fabs(fval2));
    }
    for (j = 0; j < 5; j++) {
        jtwm1 = j * 2;
        absc = hlgth * XGK21[jtwm1];
        fval1 = (*f)(centr-absc, user_data);
        fval2 = (*f)(centr+absc, user_data);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk = resk + WGK21[jtwm1] * fsum;
        *resabs = (*resabs) + WGK21[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    reskh = resk * 0.5;
    *resasc = WGK21[10] * fabs(fc - reskh);
    for (j = 0; j < 10; j++ )
        *resasc = (*resasc) + WGK21[j] * (fabs(fv1[j] - reskh) +
            fabs(fv2[j] - reskh));
    result = resk * hlgth;
    *resabs = (*resabs) * dhlgth;
    *resasc = (*resasc) * dhlgth;
    *abserr = fabs((resk - resg) * hlgth);
    if ((*resasc != 0.0) && (*abserr != 0.0))
        *abserr = (*resasc) * min(1.0,pow((200.0 * (*abserr)/(*resasc)),1.5));
    if (*resabs > uflow/(50.0 * epmach))
        *abserr = max(epmach * 50.0 * (*resabs),(*abserr));
    return result;
}




double G_K31(dq_function_type f,double a,double b,double *abserr,
    double *resabs,double *resasc, void* user_data)
{
    static double XGK31[16] = {
        0.99800229869339706029,
        0.98799251802048542849,
        0.96773907567913913426,
        0.93727339240070590431,
        0.89726453234408190088,
        0.84820658341042721620,
        0.79041850144246593297,
        0.72441773136017004742,
        0.65099674129741697053,
        0.57097217260853884754,
        0.48508186364023968069,
        0.39415134707756336990,
        0.29918000715316881217,
        0.20119409399743452230,
        0.10114206691871749903,
        0.00000000000000000000};
    static double WGK31[16] = {
        0.00537747987292334899,
        0.01500794732931612254,
        0.02546084732671532019,
        0.03534636079137584622,
        0.04458975132476487661,
        0.05348152469092808727,
        0.06200956780067064029,
        0.06985412131872825871,
        0.07684968075772037889,
        0.08308050282313302104,
        0.08856444305621177065,
        0.09312659817082532123,
        0.09664272698362367851,
        0.09917359872179195933,
        0.10076984552387559504,
        0.10133000701479154902};
    static double WG15[8] = {
        0.03075324199611726835,
        0.07036604748810812471,
        0.10715922046717193501,
        0.13957067792615431445,
        0.16626920581699393355,
        0.18616100001556221103,
        0.19843148532711157646,
        0.20257824192556127288};

    double fv1[15],fv2[15];
    double absc,centr,dhlgth;
    double fc,fsum,fval1,fval2,hlgth;
    double resg,resk,reskh,result;
    int j,jtw,jtwm1;

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    fc=(*f)(centr, user_data);
    resg = fc * WG15[7];
    resk = fc * WGK31[15];
    *resabs = fabs(resk);
    for (j = 0; j < 7; j++) {
        jtw = 2 * j + 1;
        absc = hlgth * XGK31[jtw];
        fval1 = (*f)(centr-absc, user_data);
        fval2 = (*f)(centr+absc, user_data);
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg += WG15[j] * fsum;
        resk += WGK31[jtw] * fsum;
        *resabs = *resabs + WGK31[jtw] * (fabs(fval1) + fabs(fval2));
    }
    for (j = 0; j < 8; j++) {
        jtwm1 = j * 2;
        absc = hlgth * XGK31[jtwm1];
        fval1 = (*f)(centr-absc, user_data);
        fval2 = (*f)(centr+absc, user_data);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk = resk + WGK31[jtwm1] * fsum;
        *resabs = (*resabs) + WGK31[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    reskh = resk * 0.5;
    *resasc = WGK31[15] * fabs(fc - reskh);
    for (j = 0; j < 15; j++ )
        *resasc = (*resasc) + WGK31[j] * (fabs(fv1[j] - reskh) +
            fabs(fv2[j] - reskh));
    result = resk * hlgth;
    *resabs = (*resabs) * dhlgth;
    *resasc = (*resasc) * dhlgth;
    *abserr = fabs((resk - resg) * hlgth);
    if ((*resasc != 0.0) && (*abserr != 0.0))
        *abserr = (*resasc) * min(1.0,pow((200.0 * (*abserr)/(*resasc)),1.5));
    if (*resabs > uflow/(50.0 * epmach))
        *abserr = max(epmach * 50.0 * (*resabs),(*abserr));
    return result;
}


#include <float.h>
#include <math.h>


double G_K41(dq_function_type f,double a,double b,double *abserr,
    double *resabs,double *resasc, void* user_data)
{
/* Gauss-Kronrod abscissae and weights for 41 - 20 rules */
    static double XGK41[21] = {
        0.99885903158827766384,
        0.99312859918509492479,
        0.98150787745025025919,
        0.96397192727791379127,
        0.94082263383175475352,
        0.91223442825132590587,
        0.87827681125228197608,
        0.83911697182221882339,
        0.79504142883755119835,
        0.74633190646015079261,
        0.69323765633475138481,
        0.63605368072651502545,
        0.57514044681971031534,
        0.51086700195082709800,
        0.44359317523872510320,
        0.37370608871541956067,
        0.30162786811491300432,
        0.22778585114164507808,
        0.15260546524092267551,
        0.07652652113349733375,
        0.00000000000000000000};
    static double WGK41[21] = {
        0.00307358371852053150,
        0.00860026985564294220,
        0.01462616925697125298,
        0.02038837346126652360,
        0.02588213360495115883,
        0.03128730677703279896,
        0.03660016975820079803,
        0.04166887332797368626,
        0.04643482186749767472,
        0.05094457392372869193,
        0.05519510534828599474,
        0.05911140088063957237,
        0.06265323755478116803,
        0.06583459713361842211,
        0.06864867292852161935,
        0.07105442355344406831,
        0.07303069033278666750,
        0.07458287540049918899,
        0.07570449768455667466,
        0.07637786767208073671,
        0.07660071191799965645};
    static double WG20[10] = {
        0.01761400713915211831,
        0.04060142980038694133,
        0.06267204833410906357,
        0.08327674157670474872,
        0.10193011981724043504,
        0.11819453196151841731,
        0.13168863844917662690,
        0.14209610931838205133,
        0.14917298647260374679,
        0.15275338713072585070};
    double fv1[20],fv2[20];
    double absc,centr,dhlgth;
    double fc,fsum,fval1,fval2,hlgth;
    double resg,resk,reskh,result;
    int j,jtw,jtwm1;

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    resg = 0.0;
    fc = (*f)(centr, user_data);
    resk = fc * WGK41[20];
    *resabs = fabs(resk);
    for (j = 0; j < 10; j++) {
        jtw = 2 * j + 1;
        absc = hlgth * XGK41[jtw];
        fval1 = (*f)(centr-absc, user_data);
        fval2 = (*f)(centr+absc, user_data);
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg += WG20[j] * fsum;
        resk += WGK41[jtw] * fsum;
        *resabs = *resabs + WGK41[jtw] * (fabs(fval1) + fabs(fval2));
    }
    for (j = 0; j < 10; j++) {
        jtwm1 = j * 2;
        absc = hlgth * XGK41[jtwm1];
        fval1 = (*f)(centr-absc, user_data);
        fval2 = (*f)(centr+absc, user_data);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk = resk + WGK41[jtwm1] * fsum;
        *resabs = (*resabs) + WGK41[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    reskh = resk * 0.5;
    *resasc = WGK41[20] * fabs(fc - reskh);
    for (j = 0; j < 20; j++ )
        *resasc = (*resasc) + WGK41[j] * (fabs(fv1[j] - reskh) +
            fabs(fv2[j] - reskh));
    result = resk * hlgth;
    *resabs = (*resabs) * dhlgth;
    *resasc = (*resasc) * dhlgth;
    *abserr = fabs((resk - resg) * hlgth);
    if ((*resasc != 0.0) && (*abserr != 0.0))
        *abserr = (*resasc) * min(1.0,pow((200.0 * (*abserr)/(*resasc)),1.5));
    if (*resabs > uflow/(50.0 * epmach))
        *abserr = max(epmach * 50.0 * (*resabs),(*abserr));
    return result;
}


double G_K51(dq_function_type f,double a,double b,double *abserr,
    double *resabs,double *resasc, void* user_data)
{
/* Gauss-Kronrod abscissae and weights for 51 - 25 rules */
    static double XGK51[26] = {
        0.99926210499260983419,
        0.99555696979049809791,
        0.98803579453407724764,
        0.97666392145951751150,
        0.96161498642584251242,
        0.94297457122897433941,
        0.92074711528170156175,
        0.89499199787827536885,
        0.86584706529327559545,
        0.83344262876083400142,
        0.79787379799850005941,
        0.75925926303735763058,
        0.71776640681308438819,
        0.67356636847346836449,
        0.62681009901031741279,
        0.57766293024122296772,
        0.52632528433471918260,
        0.47300273144571496052,
        0.41788538219303774885,
        0.36117230580938783774,
        0.30308953893110783017,
        0.24386688372098843205,
        0.18371893942104889202,
        0.12286469261071039639,
        0.06154448300568507889,
        0.00000000000000000000};
    static double WGK51[26] = {
        0.00198738389233031593,
        0.00556193213535671376,
        0.00947397338617415161,
        0.01323622919557167481,
        0.01684781770912829823,
        0.02043537114588283546,
        0.02400994560695321622,
        0.02747531758785173780,
        0.03079230016738748889,
        0.03400213027432933784,
        0.03711627148341554356,
        0.04008382550403238207,
        0.04287284502017004948,
        0.04550291304992178891,
        0.04798253713883671391,
        0.05027767908071567196,
        0.05236288580640747586,
        0.05425112988854549014,
        0.05595081122041231731,
        0.05743711636156783285,
        0.05868968002239420796,
        0.05972034032417405998,
        0.06053945537604586295,
        0.06112850971705304831,
        0.06147118987142531666,
        0.06158081806783293508};
    static double WG25[13] = {
        0.01139379850102628795,
        0.02635498661503213726,
        0.04093915670130631266,
        0.05490469597583519193,
        0.06803833381235691721,
        0.08014070033500101801,
        0.09102826198296364981,
        0.10053594906705064420,
        0.10851962447426365312,
        0.11485825914571164834,
        0.11945576353578477223,
        0.12224244299031004169,
        0.12317605372671545120};

    double fv1[25],fv2[25];
    double absc,centr,dhlgth;
    double fc,fsum,fval1,fval2,hlgth;
    double resg,resk,reskh,result;
    int j,jtw,jtwm1;

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    fc=(*f)(centr, user_data);
    resg = fc * WG25[12];
    resk = fc * WGK51[25];
    *resabs = fabs(resk);
    for (j = 0; j < 12; j++) {
        jtw = 2 * j + 1;
        absc = hlgth * XGK51[jtw];
        fval1 = (*f)(centr-absc, user_data);
        fval2 = (*f)(centr+absc, user_data);
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg += WG25[j] * fsum;
        resk += WGK51[jtw] * fsum;
        *resabs = *resabs + WGK51[jtw] * (fabs(fval1) + fabs(fval2));
    }
    for (j = 0; j < 13; j++) {
        jtwm1 = j * 2;
        absc = hlgth * XGK51[jtwm1];
        fval1 = (*f)(centr-absc, user_data);
        fval2 = (*f)(centr+absc, user_data);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk = resk + WGK51[jtwm1] * fsum;
        *resabs = (*resabs) + WGK51[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    reskh = resk * 0.5;
    *resasc = WGK51[25] * fabs(fc - reskh);
    for (j = 0; j < 25; j++ )
        *resasc = (*resasc) + WGK51[j] * (fabs(fv1[j] - reskh) +
            fabs(fv2[j] - reskh));
    result = resk * hlgth;
    *resabs = (*resabs) * dhlgth;
    *resasc = (*resasc) * dhlgth;
    *abserr = fabs((resk - resg) * hlgth);
    if ((*resasc != 0.0) && (*abserr != 0.0))
        *abserr = (*resasc) * min(1.0,pow((200.0 * (*abserr)/(*resasc)),1.5));
    if (*resabs > uflow/(50.0 * epmach))
        *abserr = max(epmach * 50.0 * (*resabs),(*abserr));
    return result;
}


double G_K61(dq_function_type f,double a,double b,double *abserr,
    double *resabs,double *resasc, void* user_data)
{
/* Gauss-Kronrod abscissae and weights for 61 - 30 rules */

    static double XGK61[31] = {
        0.99948441005049063757,
        0.99689348407464954027,
        0.99163099687040459486,
        0.98366812327974720997,
        0.97311632250112626837,
        0.96002186496830751222,
        0.94437444474855997942,
        0.92620004742927432588,
        0.90557330769990779855,
        0.88256053579205268154,
        0.85720523354606109896,
        0.82956576238276839744,
        0.79972783582183908301,
        0.76777743210482619492,
        0.73379006245322680473,
        0.69785049479331579693,
        0.66006106412662696137,
        0.62052618298924286114,
        0.57934523582636169176,
        0.53662414814201989926,
        0.49248046786177857499,
        0.44703376953808917678,
        0.40040125483039439254,
        0.35270472553087811347,
        0.30407320227362507737,
        0.25463692616788984644,
        0.20452511668230989144,
        0.15386991360858354696,
        0.10280693796673703015,
        0.05147184255531769583,
        0.00000000000000000000};
    static double WGK61[31] = {
        0.00138901369867700762,
        0.00389046112709988405,
        0.00663070391593129217,
        0.00927327965951776343,
        0.01182301525349634174,
        0.01436972950704580481,
        0.01692088918905327263,
        0.01941414119394238117,
        0.02182803582160919230,
        0.02419116207808060137,
        0.02650995488233310161,
        0.02875404876504129284,
        0.03090725756238776247,
        0.03298144705748372603,
        0.03497933802806002414,
        0.03688236465182122922,
        0.03867894562472759295,
        0.04037453895153595911,
        0.04196981021516424615,
        0.04345253970135606932,
        0.04481480013316266319,
        0.04605923827100698812,
        0.04718554656929915395,
        0.04818586175708712914,
        0.04905543455502977889,
        0.04979568342707420636,
        0.05040592140278234684,
        0.05088179589874960649,
        0.05122154784925877217,
        0.05142612853745902593,
        0.05149472942945156756};
    static double WG30[15] = {
        0.00796819249616660562,
        0.01846646831109095914,
        0.02878470788332336935,
        0.03879919256962704960,
        0.04840267283059405290,
        0.05749315621761906648,
        0.06597422988218049513,
        0.07375597473770520627,
        0.08075589522942021535,
        0.08689978720108297980,
        0.09212252223778612872,
        0.09636873717464425964,
        0.09959342058679526706,
        0.10176238974840550460,
        0.10285265289355884034};

    double fv1[30],fv2[30];
    double absc,centr,dhlgth;
    double fc,fsum,fval1,fval2,hlgth;
    double resg,resk,reskh,result;
    int j,jtw,jtwm1;

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    resg = 0.0;
    fc=(*f)(centr, user_data);
    resk = fc * WGK61[30];
    *resabs = fabs(resk);
    for (j = 0; j < 15; j++) {
        jtw = 2 * j + 1;
        absc = hlgth * XGK61[jtw];
        fval1 = (*f)(centr-absc, user_data);
        fval2 = (*f)(centr+absc, user_data);
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg += WG30[j] * fsum;
        resk += WGK61[jtw] * fsum;
        *resabs = *resabs + WGK61[jtw] * (fabs(fval1) + fabs(fval2));
    }
    for (j = 0; j < 15; j++) {
        jtwm1 = j * 2;
        absc = hlgth * XGK61[jtwm1];
        fval1 = (*f)(centr-absc, user_data);
        fval2 = (*f)(centr+absc, user_data);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk += WGK61[jtwm1] * fsum;
        *resabs = *resabs + WGK61[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    reskh = resk * 0.5;
    *resasc = WGK61[30] * fabs(fc - reskh);
    for (j = 0; j < 30; j++ )
        *resasc = (*resasc) + WGK61[j] * (fabs(fv1[j] - reskh) +
            fabs(fv2[j] - reskh));
    result = resk * hlgth;
    *resabs = (*resabs) * dhlgth;
    *resasc = (*resasc) * dhlgth;
    *abserr = fabs((resk - resg) * hlgth);
    if ((*resasc != 0.0) && (*abserr != 0.0))
        *abserr = (*resasc) * min(1.0,pow((200.0 * (*abserr)/(*resasc)),1.5));
    if (*resabs > uflow/(50.0 * epmach))
        *abserr = max(epmach * 50.0 * (*resabs),(*abserr));
    return result;
}


void dqmomo(double alpha,double beta,double ri[],double rj[],
    double rg[], double rh[],int wgtfunc)
{
    double alfp1,alfp2,an,anm1,betp1,betp2,ralf,rbet;
    int i,im1;

    alfp1 = alpha + 1.0;
    betp1 = beta + 1.0;
    alfp2 = alpha + 2.0;
    betp2 = beta + 2.0;
    ralf = pow(2.0,alfp1);
    rbet = pow(2.0,betp1);

/* Compute ri, rj using a forward recurrence relation. */
    ri[0] = ralf / alfp1;
    rj[0] = rbet / betp1;
    ri[1] = ri[0] * alpha / alfp2;
    rj[1] = rj[0] * beta / betp2;
    an = 2.0;
    anm1 = 1.0;
    for(i = 2;i < 25; i++) {
        ri[i] = -(ralf + an * (an - alfp2) * ri[i - 1]) /
            (anm1 * (an + alfp1));
        rj[i] = -(rbet + an * (an - betp2) * rj[i - 1]) /
            (anm1 * (an + betp1));
        anm1 = an;
        an += 1.0;
    }
    if (wgtfunc == 1)
        goto _70;
    if (wgtfunc == 3)
        goto _40;

/* Compute rg using a forward recurrence formula. */
    rg[0] = -ri[0] / alfp1;
    rg[1] = -(ralf + ralf) / (alfp2 * alfp2) - rg[0];
    an = 2.0;
    im1 = 1;    /* FORTRAN uses im1 = 2 */
    for (i = 2; i < 25; i++) {
        rg[i] = -(an * (an - alfp2) * rg[im1] - an * ri[im1] +
            anm1 * ri[i]) / (anm1 * (an + alfp1));
        anm1 = an;
        an += 1.0;
        im1 = i;
    }
    if (wgtfunc == 2)
        goto _70;

/* Compute rh using a forward recurrence relation. */
_40:
    rh[0] = -rj[0] / betp1;
    rh[1] = -(rbet + rbet) / (betp2 * betp2) - rh[0];
    an = 2.0;
    anm1 = 1.0;
    im1 = 1;    /* FORTRAN uses im1 = 2 */
    for (i = 2; i < 25; i++) {
        rj[i] = -(an * (an - betp2) * rh[im1] - an * rj    [im1] +
            anm1 * rj[i]) / ( anm1 * (an + betp1));
        anm1 = an;
        an += 1.0;
        im1 = i;
    }
    for (i = 1; i < 25; i += 2) {
        rh[i] = -rh[i];
    }
_70:
    for (i = 1; i < 25; i += 2) {
        rj[i] = -rj[i];
    }
}




double dqng(dq_function_type f,double a,double b,double epsabs,
    double epsrel,double *abserr,int *neval,int *ier, void* user_data)
{
    static double x1[5] = {
        0.97390652851717172008,
        0.86506336668898451073,
        0.67940956829902440623,
        0.43339539412924719080,
        0.14887433898163121088};
    static double w10[5] = {
        0.06667134430868813759,
        0.14945134915058059315,
        0.21908636251598204400,
        0.26926671930999635509,
        0.29552422471475287017};
    static double x2[5] ={
        0.99565716302580808074,
        0.93015749135570822600,
        0.78081772658641689706,
        0.56275713466860468334,
        0.29439286270146019813};
    static double w21a[5] = {
        0.03255816230796472748,
        0.07503967481091995277,
        0.10938715880229764190,
        0.13470921731147332593,
        0.14773910490133849137};
    static double w21b[6] = {
        0.01169463886737187428,
        0.05475589657435199603,
        0.09312545458369760554,
        0.12349197626206585108,
        0.14277593857706008080,
        0.14944555400291690566};
    static double x3[11] = {
        0.99933336090193208139,
        0.98743340290808886980,
        0.95480793481426629926,
        0.90014869574832829363,
        0.82519831498311415085,
        0.73214838898930498261,
        0.62284797053772523864,
        0.49947957407105649995,
        0.36490166134658076804,
        0.22225491977660129650,
        0.07465061746138332204};
    static double w43a[10] = {
        0.01629673428966656492,
        0.03752287612086950146,
        0.05469490205825544215,
        0.06735541460947808608,
        0.07387019963239395343,
        0.00576855605976979618,
        0.02737189059324884208,
        0.04656082691042883074,
        0.06174499520144256450,
        0.07138726726869339777};
    static double w43b[12] = {
        0.00184447764021241410,
        0.01079868958589165174,
        0.02189536386779542810,
        0.03259746397534568944,
        0.04216313793519181185,
        0.05074193960018457778,
        0.05837939554261924838,
        0.06474640495144588554,
        0.06956619791235648453,
        0.07282444147183320815,
        0.07450775101417511827,
        0.07472214751740300559};
    static double x4[22] = {
        0.99990297726272923449,
        0.99798989598667874543,
        0.99217549786068722281,
        0.98135816357271277357,
        0.96505762385838461913,
        0.94316761313367059682,
        0.91580641468550720959,
        0.88322165777131650137,
        0.84571074846241566661,
        0.80355765803523098279,
        0.75700573068549555833,
        0.70627320978732181982,
        0.65158946650117792253,
        0.59322337405796108888,
        0.53149360597083193229,
        0.46676362304202284487,
        0.39942484785921880473,
        0.32987487710618828827,
        0.25850355920216155180,
        0.18569539656834665202,
        0.11184221317990746817,
        0.03735212339461987081};
    static double w87a[21] = {
        0.00814837738414917290,
        0.01876143820156282224,
        0.02734745105005228616,
        0.03367770731163793005,
        0.03693509982042790761,
        0.00288487243021153050,
        0.01368594602271270189,
        0.02328041350288831112,
        0.03087249761171335868,
        0.03569363363941877072,
        0.00091528334520224136,
        0.00539928021930047137,
        0.01094767960111893113,
        0.01629873169678733526,
        0.02108156888920383511,
        0.02537096976925382724,
        0.02918969775647575250,
        0.03237320246720278969,
        0.03478309895036514275,
        0.03641222073135178756,
        0.03725387550304770854};
    static double w87b[23] = {
        0.00027414556376207235,
        0.00180712415505794295,
        0.00409686928275916486,
        0.00675829005184737870,
        0.00954995767220164654,
        0.01232944765224485369,
        0.01501044734638895238,
        0.01754896798624319110,
        0.01993803778644088820,
        0.02219493596101228680,
        0.02433914712600080547,
        0.02637450541483920724,
        0.02828691078877120066,
        0.03005258112809269532,
        0.03164675137143992940,
        0.03305041341997850329,
        0.03425509970422606179,
        0.03526241266015668103,
        0.03607698962288870119,
        0.03669860449845609450,
        0.03712054926983257611,
        0.03733422875193504032,
        0.03736107376267902341};
    double fv1[5],fv2[5],fv3[5],fv4[5],savfun[21];
    double absc,centr,dhlgth;
    double fcentr,fval,fval1,fval2,hlgth;
    double result,res10,res21,res43,res87;
    double resabs,resasc,reskh;
    int ipx,k,l;

    result = 0.0;
    *abserr = 0.0;
    *neval = 0;
    *ier = 6;
    if ((epsabs < 0.0) && (epsrel < 0.0)) return result;
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);
    centr = 0.5 * (a + b);
    fcentr=(*f)(centr, user_data);
    *neval = 21;
    *ier = 1;

    for (l = 1; l <=3; l++) {
        switch (l) {
            case 1:
                res10 = 0.0;
                res21 = w21b[5] * fcentr;
                resabs = w21b[5] * fabs(fcentr);
                for (k = 0;k < 5; k++) {
                    absc = hlgth * x1[k];
                    fval1 = (*f)(centr+absc, user_data);
                    fval2 = (*f)(centr-absc, user_data);
                    fval = fval1 + fval2;
                    res10 += (w10[k] * fval);
                    res21 += (w21a[k] * fval);
                    resabs += (w21a[k] *
                        (fabs(fval1) + fabs(fval2)));
                    savfun[k] = fval;
                    fv1[k] = fval1;
                    fv2[k] = fval2;
                }
                ipx = 4;
                for (k = 0; k < 5; k++) {
                    ipx++;
                    absc = hlgth * x2[k];
                    fval1 = (*f)(centr + absc, user_data);
                    fval2 = (*f)(centr - absc, user_data);
                    fval = fval1 + fval2;
                    res21 += (w21b[k] * fval);
                    resabs += (w21b[k] *
                        (fabs(fval1) + fabs(fval2)));
                    savfun[ipx] = fval;
                    fv3[k] = fval1;
                    fv4[k] = fval2;
                }
                result = res21 * hlgth;
                resabs *= dhlgth;
                reskh = 0.5 * res21;
                resasc = w21b[5] * fabs(fcentr - reskh);
                for (k = 0; k < 5; k++)
                 resasc += (w21a[k] * (fabs(fv1[k] -reskh) +
                  fabs(fv2[k] - reskh)) + w21b[k] *
                  (fabs(fv3[k] - reskh) + fabs(fv4[k]-reskh)));
                *abserr = fabs((res21 - res10) * hlgth);
                resasc *=dhlgth;
                break;
            case 2:
                res43 = w43b[11] * fcentr;
                *neval = 43;
                for (k = 0; k < 10; k++)
                 res43 += (savfun[k] * w43a[k]);
                for (k = 0; k < 11; k++) {
                    ipx++;
                    absc = hlgth * x3[k];
                    fval = (*f)(centr+absc, user_data) + (*f)(centr-absc, user_data);
                    res43 += (fval * w43b[k]);
                    savfun[ipx] = fval;
                }
                result = res43 * hlgth;
                *abserr = fabs((res43-res21) * hlgth);
                break;
            case 3:
                res87 = w87b[22] * fcentr;
                *neval = 87;
                for (k = 0; k < 21; k++)
                    res87 += (savfun[k] * w87a[k]);
                for (k = 0; k < 22; k++) {
                    absc = hlgth * x4[k];
                    res87 += w87b[k] *
                        ((*f)(centr+absc, user_data) + (*f)(centr-absc, user_data));
                }
                result = res87 * hlgth;
                *abserr = fabs((res87 - res43) * hlgth);
        }
        if ((resasc != 0.0) &&(*abserr != 0.0))
            *abserr = resasc * min(1.0,
                pow(200.0 * (*abserr)/resasc,1.5));
        if (resabs > uflow/(50.0 * epmach)) *abserr =
            max((epmach * 50.0) * resabs,*abserr);
        if (*abserr <= max(epsabs,epsrel * fabs(result))) *ier = 0;
        if (*ier == 0) break;
    }
    return result;
}


void dqsort(int limit,int last,int *maxerr,double *ermax,double elist[],
    int iord[],int *nrmax)
{
    double errmax,errmin;
    int i,ibeg,ido,isucc,j,jbnd,jupbn,k;

    if (last > 1) goto _10;
    iord[0] = 0;
    iord[1] = 1;
    goto _90;
_10:
    errmax = elist[*maxerr];
    if (*nrmax == 0) goto _30;
    ido = (*nrmax) - 1;
    for (i = 0;i <= ido; i++) {
        isucc = iord[*nrmax-1];
        if (errmax <= elist[isucc]) goto _30;
        iord[*nrmax] = isucc;
        (*nrmax)--;
    }
_30:
    jupbn = last;
    if (last > (limit/2 + 2))
        jupbn = limit + 3 - last;
    errmin = elist[last];
    jbnd = jupbn - 1;
    ibeg = *nrmax + 1;
    if (ibeg > jbnd) goto _50;
    for (i = ibeg; i <= jbnd; i++) {
        isucc = iord[i];
        if (errmax >= elist[isucc]) goto _60;
        iord[i-1] = isucc;
    }
_50:
    iord[jbnd] = *maxerr;
    iord[jupbn] = last;
    goto _90;
_60:
    iord[i-1] = *maxerr;
    k = jbnd;
    for (j = i;j <= jbnd; j++) {
        isucc = iord[k];
        if (errmin < elist[isucc]) goto _80;
        iord[k+1] = isucc;
        k--;
    }
    iord[i] = last;
    goto _90;
_80:
    iord[k+1] = last;
_90:
    *maxerr = iord[*nrmax];
    *ermax = elist[*maxerr];
    return;
}


double dqwgtc(double x,double c,double p2,double p3,double p4,int kp)
{
    return 1.0/(x - c);
}

double dqwgto(double x,double omega,double p2,double p3,double p4,int wgtfunc)
{
    double omx;

    omx = omega * x;
    if (wgtfunc == 1)
        return cos(omx);
    else
        return sin(omx);
}

double dqwgts(double x,double a,double b,double alpha,double beta,int wgtfunc)
{
    double bmx,xma,result;

    xma = x - a;
    bmx = b - x;
    result = pow(xma,alpha) * pow(bmx,beta);
    switch (wgtfunc) {
        case 1:
            return result;
        case 2:
            return result * log(xma);
        case 3:
            return result * log(bmx);
        case 4:
            return result * log(xma) * log(bmx);
        default:
            return result;
    }
}


/*
          EXAMPLE of 2D Integration
*/


/*

       dblint.c -- example of double integrator.
 *
 *      Main routine calls the outer routine which points
 *      to a function which computes the upper and lower
 *		integration limits for the inner routine and
 *      returns the computed integral.
 
 		
#include <stdio.h>
#include <float.h>
#include <math.h>
#include "cquadpack.h"

double y_global;

double f(double x)
{
	return log(x * x * y_global * y_global);
}

double f1(double x)
{
    return x * x * y_global;
}

double g(double x)
{
	double y,a,b,abserr,resabs,resasc;

	a = 1.0;
	b = sqrt(2.0);
	y_global = x;
	y = G_K61(f,a,b,&abserr,&resabs,&resasc);
	return y;
}
double g1(double x)
{
    double y,a,b,abserr,resabs,resasc;

    a = 0.0;
    b = 1.0;
    y_global = x;   
    y = G_K61(f1,a,b,&abserr,&resabs,&resasc);
    return y;
}

void main()
{
    double y,a1,b1,abserr,resabs,resasc;
	
    a1 = 0;
    b1 = 1;

    y = G_K61(g1,a1,b1,&abserr,&resabs,&resasc);
	
	printf("Integral = %.15lg\n",y);
	printf("abserr = %.15\lg\n",abserr);
}

!====================================================================================!
!             SECOND EXAMPLE
!====================================================================================!
   

    #include <stdio.h>
#include "cquadpack.h"

double x_global;

double g(x)
{
    return x;
}
double h(x)
{
    return x*x;
}

double z(double f(),double x,double g(),double h(),
    double *abserr,double *resabs,double *resasc)
{
    double c,d,result;
    c = g(x);
    d = h(x);
    x_global = x;
    result = G_K61(f,c,d,abserr,resabs,resasc);

}
double f(double y)
{
    return x_global * x_global * y;
}

double dblint(double f(),double a,double b,double g(),double h(),
    double *abserr,double *resabs,double *resasc)
{
    double result;

    result = G_K61(z,a,b,abserr,resabs,resask);
   
    
}
void main()
{
    double a,b;

    a = 0.0;
    b = 1.0;
    result = dblint(f,a,b,g,h,&abserr);
    printf ("Result = %.12lg\n",result);
}

!====================================================================!
!                  THIRD EXAMPLE
!====================================================================!
      Using quadpack to compute elliptic integral.
#include <stdio.h>
#include <float.h>
#include <math.h>
#include "cquadpack.h"

double EllipticE(double x)
{
	return sqrt(1.0-0.5*sin(x)*sin(x));
}

void main()
{
	double y,a,b,abserr,resabs,resasc;
	
	a = 0.0;
	b = M_PI_2;
	y = G_K21(EllipticE,a,b,&abserr,&resabs,&resasc);
	printf("Integral = %.15lg\n",y);
	printf("abserr = %.15\lg\n",abserr);
}

*/







