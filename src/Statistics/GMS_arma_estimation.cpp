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




#if defined __GNUC__ && !defined __INTEL_COMPILER
   #include <omp.h>
#endif
#include <algorithm>
#include "GMS_arma_estimation.h"
#include "GMS_cephes.h"

/*
      Converted to C++ from F90 implementation of the following ASA algorithms
      ALGORITHM AS 154  APPL. STATIST. (1980) VOL.29, P.311
      ALGORITHM AS 182  APPL. STATIST. (1982) VOL.31, NO.2
      Modified by Bernard Gingold on 29-11-2020 SUN 2:53PM +00200
*/

/*
   !  ALGORITHM AS 182  APPL. STATIST. (1982) VOL.31, NO.2

!  Finite sample prediction from ARIMA processes.
*/



void gms::stat::forkal(const int32_t ip,
            const int32_t iq,
            const int32_t ir,
            const int32_t np,
            const int32_t ird,
            const int32_t irz,
            const int32_t id,
            const int32_t il,
            const int32_t n,
            const int32_t nrbar,
            float * __restrict __attribute__((aligned(64))) phi,
            float * __restrict __attribute__((aligned(64))) theta,
            float * __restrict __attribute__((aligned(64))) delta,
            float * __restrict __attribute__((aligned(64))) w,
            float * __restrict __attribute__((aligned(64))) y,
            float * __restrict __attribute__((aligned(64))) amse,
            float * __restrict __attribute__((aligned(64))) a,
            float * __restrict __attribute__((aligned(64))) p,
            float * __restrict __attribute__((aligned(64))) v,
            float * __restrict __attribute__((aligned(64))) resid,
            float * __restrict __attribute__((aligned(64))) e,
            float * __restrict __attribute__((aligned(64))) xnext,
            float * __restrict __attribute__((aligned(64))) xrow,
            float * __restrict __attribute__((aligned(64))) rbar,
            float * __restrict __attribute__((aligned(64))) thetab,
            float * __restrict __attribute__((aligned(64))) store,
            int32_t &ifault ) {
      //
    // Invoking this routine will calculate the finite sample predictions
     //     and their conditional mean square errors for any ARIMA process.
     // Will surely impact the branch predictor resources.
     int32_t k;
     ifault = 0;
     if(ip<0) ifault = 1;
     if(iq<0) ifault = ifault + 2;
     if(ip*ip+iq*iq == 0) ifault = 4;
     k = iq + 1;
     if(k<ip) k = ip;
     if(ir!=k) ifault = 5;
     if(np!=ir*(ir+1)/2) ifault = 6;
     if(nrbar!=np*(np-1)/2) ifault = 7;
     if(id<0) ifault = 8;
     if(ird!=ir+id) ifault = 9;
     if(irz!=ird*(ird+1)/2) ifault = 10;
     if(il<1) ifault = 11;
     if(ifault!=0) return;

     int32_t i, i45, ibc, id1, id2r, id2r1, id2r2, idd1, idd2, iddr, idrr1, 
            idk, iid, ind, ind1, ind2, iq1, ir1, ir2, iri, iri1, irj, iupd, 
            j, j1, jj, jkl, jkl1, jklj, jrj, jrk, k, k1, kk, kk1, kkk, l,
            lk, lk1, ll, lli, nit, nj, nt;
     float a1, aa, ams, del, dt, phii, phij, phijdt, sigma, ssq, sumlog;
     //  Calculate initial conditions for Kalman filter
#if defined __INTEL_COMPILER || defined __ICC
     __assume_aligned(a,64);
     __assume_aligned(v,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
     a = (float*)__builtin_assume_aligned(a,64);
     v = (float*)__builtin_assume_aligned(v,64);
#endif
     a[0] = 0.0f;
     v[0] = 0.0f;
     if(np==1) goto L130;
#if defined __INTEL_COMPILER
#pragma vector aligned
#pragma vector always
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(v:64)
#endif
     for(int32_t idx = 2; idx != np; ++idx) {
         v[i] = 0.0f;
     }
     if(iq==0) goto L130;
     iq1 = iq+1;
#if defined __INTEL_COMPILER || defined __ICC
     __assume_aligned(theta,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
     theta = (float*)__builtine_assume_aligned(theta,64);
#endif
#if defined __INTEL_COMPILER
#pragma ivdep
#pragma vector aligned
#pragma vector always
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(theta:64,v:64)
#endif
     for(i = 1; i != iq1; ++i) {
         v[i] = theta[i-1];
     }
     for(j = 0; j != iq; ++j) { // Not easy vectorizable (will probably generate a gether-scatter stor/loads for index lli)
        ll = j * (ir+ir+1-j) / 2;
        float tj = theta[j];
        for(i = j; i != iq; ++i) {
            lli = ll + i;
            v[lli] = theta[i] * tj;
        }
    }
    /*
         Find initial likelihood conditions.
!     IFAULT not tested on exit from STARMA as all possible errors
!     have been checked above.
    */
L130:
     if(ir==1) p[1] = 1.0f/(1.0f-phi[0]*phi[0]);
     if(ir!=1) {
         starma(ip,iq,ir,np,phi,theta,a,p,v,
                thetab,xnext,xrow,rbar,nrbar,ifault);
     }  
     //    Calculate data transformations
     nt = n - id;
     if(id==0) goto L170;
#if defined __INTEL_COMPILER || defined __ICC
     __assume_aligned(store,64);
     __assume_aligned(w,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
     store = (float*)__builtin_assume_aligned(store,64);
     w     = (float*)__builtin_assume_aligned(w,64);
#endif
#if defined __INTEL_COMPILER || defined __ICC
#pragma vector aligned
#pragma vector always
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(store:64,w:64)
#endif
     for(j = 1; j != id; ++j) {
         nj = n - j;
         store[j] = w[nj];
     }
#if defined __INTEL_COMPILER || defined __ICC
     __assume_aligned(delta,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
     delta = (float*)__builtin_assume_aligned(delta,64);
#endif
     for(i = 0; i != nt; ++i) {
         aa = 0.0f;
         for(k = 0; k != id; ++k) {
             idk = id + i - k;
             aa = aa - delta[k] * w[idk];
         }
         iid = i + id;
         w[i] = w[iid] + aa;
     }
     //  Evaluate likelihood to obtain final KF conditions
L170:
     sumlog = 0.0f;
     ssq    = 0.0f;
     iupd   = 1;
     del    = -1.0f;
     nit    = 0;
     karma(ip,iq,ir,phi,theta,a,p,v,nt,w,resid,
           sumlog,ssq,iupd,del,e,nit);
     //  Calculate M.L.E. of sigma squared
     sigma = 0.0f;
#if defined __INTEL_COMPILER || defined __ICC
     __assume_aligned(resid,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
    resid = (float*)__builtin_assume_aligned(resid,64);
#endif
#if defined __INTEL_COMPILER || defined __ICC
#pragma vector aligned
#pragma simd reduction(+:sigma)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd reduction(+:sigma) aligned(resid:64)
#endif
     for(j = 0; j != nt; ++j) {
         sigma = sigma+resid[j]*resid[j];
     }
     // Reset the initial A and P when differencing occurs
     if(id==0) goto L250;
#if defined __INTEL_COMPILER || defined __ICC
     __assume_aligned(xrow,64);
     __assume_aligned(p,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
     xrow = (float*)__builtin_assume_aligned(xrow,64);
     p    = (float*)__builtin_assume_aligned(p,64);
#endif
#if defined __INTEL_COMPILER || defined __ICC
#pragma vector aligned
#pragma vector always
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(xrow:64,p:64)
#endif
     for(i = 0; i != np; ++i) {
         xrow[i] = p[i];
     }
#if defined __INTEL_COMPILER || defined __ICC
#pragma vector aligned
#pragma vector always
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(p:64)
#endif
     for(i = 0; i != irz; ++i) {
         p[i] = 0.0f;
     } 
     ind = 0;
     // Probably not easily vectorizable (leave it as a scalar loop)
     for(j = 1; j != ir; ++j) {
         k = (j-1) * (id + ir + 1) - (j-1) * j / 2;
	 for(i = j; i != ir; ++i) {
             ind += 1;
	     k   += 1;
	     p[k] = xrow[ind];
	 }
     }
#if defined __INTEL_COMPILER || defined __ICC
#pragma vector aligned
#pragma vector always
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(a:64,store:64)
#endif
     for(j = 0; j != id; ++j) {
         irj = ir + 1;
	 a[irj] = store[j];
     }
     // Set up constants
L250:
     ir2 = ir + 1;
     ir1 = ir - 1;
     id1 = id - 1;
     id2r = 2 * ird;
     id2r1 = id2r - 1;
     idd1 = 2 * id + 1;
     idd2 = idd1 + 1;
     i45 = id2r + 1;
     idrr1 = ird + 1;
     iddr = 2 * id + ir;
     jkl = ir * (iddr + 1) / 2;
     jkl1 = jkl + 1;
     id2r2 = id2r + 2;
     ibc = ir * (i45 - ir) / 2;
#if defined __INTEL_COMPILER || defined __ICC
     assume_aligned(asme,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
     asme = (float*)__builtin_assume_aligned(asme,64);
#endif
     // Large processing loop
     for(l = 0; l != il; ++l) {
     // Predict A
#if defined __INTEL_COMPILER || defined __ICC
     __assume_aligned(phi,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
     phi = (float*)__builtin_assume_aligned(phi,64);
#endif
     a1 = a[0];
     if(ir==1) goto L310;
     for(i = 0; i != ir1; ++i) {
        register float t = a[i+1];
	a[i] = t;
     }
L310:
     a[ir-1] = 0.0f;
     if(ip==0) goto L330;
#if defined __INTEL_COMPILER || defined __ICC
#pragma vector aligned
#pragma vector always
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(a:64,phi:64)
#endif
      for(j = 0; j != ip; ++j) {
          a[j] = a[j]+phi[j]*a1;
      }
L330:
     if(id==0) goto L360;
#if defined __INTEL_COMPILER || defined __ICC
#pragma vector aligned
#pragma simd reduction(+:a1)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd reduction(+:a1) aligned(a:64,delta:64)
#endif
      for(j = 1; j != id; ++j) {
          irj = ir+j;
	  a1 = a1+delta[j]*a[irj];
      }
      if(id<2) goto L360;
      for(i = 1; i != id1; ++i) {
          iri1 = ird-i;
	  a[iri1+1] = a[iri1];
      }
L360:
      a[ir2] = a1;
      // Predict P
      if(id==0) goto L480;
      for(i = 0; i != id; ++i) {
         store[i] = 0.0f;
	 for(j = 0; j != id; ++j) {
             ll = std::max(i,j);
	     k  = std::min(i,j);
	     jj = jkl + (ll - k) + 1 + (k-1) * (idd2 - k) / 2;
	     store[i] = store[i] + delta[j] * p[jj];
	 }
      }

      for(j = 0; j != id1; ++j) {
          jj = id - j;
	  lk = (jj-1) * (idd2 - jj) / 2 + jkl;
          lk1 = jj * (idd1 - jj) / 2 + jkl;
	  for(i = 0; i != j; ++i) {
              lk  += 1;
	      lk1 += 1;
	      p[lk1] = p[lk];
	  }
      }
      for(j = 0; j != id1; ++j) {
          jklj = jkl1 + j;
	  irj  = ir + j;
	  p[jklj] = store[j] + p[irj];
      }
      p[jkl1] = p[0];
      for(i = 0; i != id; ++i) {
          iri = ir + i;
	  p[jkl1] = p[jkl1] + delta[i] * (store[i] + 2.0f * p[iri]);
      }
#if defined __INTEL_COMPILER || defined __ICC
#pragma vector aligned
#pragma vector always
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(store:64,p:64)
#endif
     for(i = 0; i != id; ++i) {
         iri = ir + i;
	 store[i] = p[iri];
     }
     for(j = 1; j != ir; ++j) {
          kk1 = j * (id2r1 - j) / 2 + ir;
          k1 = (j-1) * (id2r - j) / 2 + ir;
	  for(i = 0; i != id; ++i) {
                kk = kk1 + i;
                k = k1 + i;
		p[k] = phi[j] * store[i];
		if(j!=ir) p[k] = p[k]+p[kk];
	  }
     }
     for(j = 1; j != ir; ++j) {
         store[j] = 0.0f;
	 kkk = j * (i45 - j) / 2 - id;
	 for( i = 0; i != id; ++i) {
              kkk += 1;
	      store[j] = store[j] + delta[i] * p[kkk];
	 }
     }
     if(id==1) goto L460;
     for(j = 1; j != ir; ++j) {
          k = j * idrr1 - j * (j+1) / 2 + 1;
	  for(i = 1; i != id1; ++i) {
              k -= 1;
	      p[k] = p[k-1];
	  }
     }
L460:
     for(j = 1; j != ir; ++j) {
         k = (j-1) * (id2r - j) / 2 + ir + 1;
	 p[k] = store[j]+phi[j]*p[0];
	 if(j<ir) p[k] += p[j+1];
     }
L480:
#if defined __INTEL_COMPILER || defined __ICC
#pragma vector aligned
#pragma vector always
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(store:64,p:64)
#endif
        for(int32_t i = 0; i != ir; ++i) {
            store[i] = p[i];
	}
	ind = 0;
	dt = p[0];
	for(j = 1; j != ir; ++j) {
            phij = phi[j];
	    phijdt = phij * dt;
	    ind2 = (j-1) * (id2r2 - j) / 2;
            ind1 = j * (i45 - j) / 2;
	    for(i = j; i != ir; ++i) {
                ind += 1;
		ind2 += 1;
		phii = phi[i];
		p[ind2] = v[ind]+phii*phijdt;
		if(j<ir) p[ind2] = p[ind2]+store[j+1]*phii;
		if(i==ir) continue;
		ind1 += 1;
		p[ind2] = p[ind2]+store[i+1]*phij+p[ind1];
	    }
	}
	// Predict Y
#if defined __INTEL_COMPILER || defined __ICC
        __assume_aligned(y,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        y = (float*)__builtin_assume_aligned(y,64)
#endif
        y[l] = a[0];
        if(id==0) goto L520:
#if defined __INTEL_COMPILER || defined __ICC
#pragma vector aligned
#pragma vector always
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(store:64,p:64)
#endif
        for(j = 0; j != id; ++j) {
            irj = ir+j;
	    y[l] = y[l] + a[irj] * delta[j];
	}
	//  Calculate M.S.E. of Y
L520:
        ams = p[0];
	for(j = 1; j != id; ++j) {
             jrj = ibc + (j-1) * (idd2 - j) / 2;
             irj = ir + j;
	     ams = ams + 2.0f * delta[j] * p[irj] + p[jrj+1] * delta[j] * delta[j];
	}
	for(j = 1; j != id1; ++j) {
            j1 += 1;
	    jrk = ibc + 1 + (j-1) * (idd2 - j) / 2;
	    for(i = j1; i != id; ++i) {
                jrk += 1;
		ams = ams + 2.0f * delta[i] * delta[j] * p[jrk];
	    }
	}
	amse[l] = ams * sigma;
     }
}

/*
    !  INVOKING THIS SUBROUTINE SETS THE VALUES OF V AND PHI, AND
!  OBTAINS THE INITIAL VALUES OF A AND P.
!  THIS ROUTINE IS NOT SUITABLE FOR USE WITH AN AR(1) PROCESS.
!  IN THIS CASE THE FOLLOWING INSTRUCTIONS SHOULD BE USED FOR INITIALISATION.
!  V(1) = 1.0
!  A(1) = 0.0
!  P(1) = 1.0 / (1.0 - PHI(1) * PHI(1))
*/


void gms::stat::starma(const int32_t ip,
            const int32_t iq,
            const int32_t ir,
            const int32_t np,
            float * __restrict __attribute__((aligned(64))) phi,
            float * __restrict __attribute__((aligned(64))) theta,
            float * __restrict __attribute__((aligned(64))) a,
            float * __restrict __attribute__((aligned(64))) p,
            float * __restrict __attribute__((aligned(64))) v,
            float * __restrict __attribute__((aligned(64))) thetab,
            float * __restrict __attribute__((aligned(64))) xnext,
            float * __restrict __attribute__((aligned(64))) xrow,
            float * __restrict __attribute__((aligned(64))) rbar,
            const int32_t nrbar,
            int32_t &ifault) {
    //  CHECK FOR FAILURE INDICATION.
    // The chain of 10 branches likely a performance hit.
    // THere should not be more than 15 branches in chain of
    // instructions.
    // Forward branches predicted not-taken.
    int32_t k;
    ifault = 0;
    if(ip<0) { ifault = 1;} // likely not-taken
    if(iq<0) { ifault += 2;} // likely not-taken
    if(ip==0 && iq==0) { ifault = 4;} // likely not-taken
    k = iq + 1;
    if(k<ip) { k = ip;}
    if(ir!=k) { ifault = 5;}
    if(np!=ir*(ir+1)/2) { ifault = 6;}
    if(nrbar!=np*(np-1)/2) { ifault = 7;}
    if(ir==1) { ifault = 8;}
    if(ifault!=0) { return;}
    //loat vj,phii,phij,ssqerr,recres,ynext;
    register float phii,phij,ynext;
    float ssqer,recres;
    int32_t ifail,irank,npr,npr1
    register int32_t ind,ind1,ind2,indj,indi;
    // NOW SET A(0), V AND PHI.
#if defined __INTEL_COMPILER || defined __ICC
    __assume_aligned(phi,64);
    __assume_aligned(theta,64);
    __assume_aligned(a,64);
    __assume_aligned(p,64);
    __assume_aligned(v,64);
    __assume_aligned(thetab,64);
    __assume_aligned(xnext,64);
    __assume_aligned(xrow,64);
    __assume_aligned(rbar,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
    phi    = (float*)__builtin_assume_aligned(phi,64);
    theta  = (float*)__builtin_assume_aligned(theta,64);
    a      = (float*)__builtin_assume_aligned(a,64);
    p      = (float*)__builtin_assume_aligned(p,64);
    v      = (float*)__builtin_assume_aligned(v,64);
    thetab = (float*)__builtin_assume_aligned(thetab,64);
    xnext  = (float*)__builtin_assume_aligned(xnext,64);
    xrow   = (float*)__builtin_assume_aligned(xrow,64);
    rbar   = (float*)__builtin_assume_aligned(rbar,64);
#endif
#if defined __INTEL_COMPILER || defined __ICC
#pragma ivdep
#pragma vector aligned
#pragma vector always
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(a:64,v:64,theta:64)
#endif
     for(int32_t i = 1; i != ir; ++i) {
         a[i] = 0.0f;
         if(i>ip) phi[i] = 0.0f;
         v[i] = 0.0f;
         if(i<=iq+1) v[i] = theta[i-1];
     }   
     a[0] = 0.0f;
     if(ip==0) phi[0] = 0.0f;
     v[0] = 0.0f;
     ind = ir;
     for(int32_t j = 1; j != ir; ++j) {
         vj = v[j];
//#if defined __INTEL_COMPILER || defined __ICC
//#pragma vector aligned
//#pragma vector always
//#elif defined __GNUC__ && !defined __INTEL_COMPILER
//#pragma omp simd aligned(v:64)
//#endif
        for(int32_t i = j; i != ir; ++i) {
            ind += 1;
            v[ind] = v[i] * vj;
        }
     }

    if(ip!=0) {
     //  NOW FIND P(0).
     /*
    
         THE SET OF EQUATIONS S * VEC(P(0)) = VEC(V) IS SOLVED FOR VEC(P(0)).
!   S IS GENERATED ROW BY ROW IN THE ARRAY XNEXT.
!   THE ORDER OF ELEMENTS IN P IS CHANGED, SO AS TO
!   BRING MORE LEADING ZEROS INTO THE ROWS OF S,
!   HENCE ACHIEVING A REDUCTION OF COMPUTING TIME.
     */
          irank = 0;
          ssqer = 0.0f;
#if defined __INTEL_COMPILER || defined __ICC
#pragma vector aligned
#pragma vector always
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(rbar:64)
#endif
        for(int32_t i = 0; i != nrbar; ++i) { rbar[i] = 0.0f;}
#if defined __INTEL_COMPILER || defined __ICC	
#pragma vector aligned
#pragma vector always
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(p:64,thetab:64,xnext:64)
#endif
            for(int32_t i = 0; i != np; ++i) {
                p[i] = 0.0f;
                thetab[i] = 0.0f;
                xnext[i] = 0.0f;
            }
            ind  = 0;
            ind1 = 0;
            npr  = np - ir;
            npr1 = npr + 1;
            indj = npr1;
            ind2 = npr;
            for(int32_t j = 0; j != ir; ++j) {
                phij = phi[j];
                xnext[indj] = 0.0f;
                indj += 1;
                indi = npr1 + j;
                for(int32_t i = j; i != ir; ++i) {
                     ind += 1;
                     ynext = v[ind];
                     phii = phi[i];
                     if(j != ir) {
                        xnext[indj] = -phii;
                        if(i != ir) {
                           xnexr[indi] = xnext[indi] - phij;
                           ind1 += 1;
                           xnext[ind1] = -1.0f;
                        }
                 }
                 xnext[npr1] = -phii * phij;
                 ind2 = ind2 + 1;
                 if(ind2>np) ind2 = 1;
                 xnext[ind2] = xnext[ind2] + 1.0f;
                 inclu2(np,1.0f,xnext,xrow,ynext,p,rbar,thetab,
                   ssqer,recres,irank,ifault);
                 // NO NEED TO CHECK IFAIL AS WEIGHT = 1.0
                 xnext[ind2] = 0.0f;
                 if(i != ir) {
                    xnext[indi] = 0.0f;
                    indi += 1;
                    xnext[ind1] = 0.0f;
                }
            }
        }
          regres(np,nrbar,rbar,thetab,p);
        // NOW RE-ORDER P.
          ind = npr;
#if defined __INTEL_COMPILER || defined __ICC	
#pragma vector aligned
#pragma vector always
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(p:64,xnext:64)
#endif
          for(int32_t i = 0; i != ir; ++i) {
              ind += 1;
              xnext[i] = p[ind];
           }
           ind = np;
           ind1 = npr;
           for(int3_t i = 0; i != npr; ++i) {
               p[ind] = p[ind1];
               ind  -= 1;
               ind1 -= 1;
           }
#if defined __INTEL_COMPILER || defined __ICC	
#pragma vector aligned
#pragma vector always
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(p:64,xnext:64)
#endif
           for(int32_t i = 0; i != ir; ++i) {
               p[i] = xnext[i];
           } 
           return;
       }
       //  P(0) IS OBTAINED BY BACKSUBSTITUTION FOR A MOVING AVERAGE PROCESS.
       indn = np+1;
       ind  = np+1;
       for(int32_t i = 0; i != ir; ++i) {
           for(int32_t j = 0; j != i; ++j) {
               ind -= 1;
               p[ind] = v[ind];
               if(j != 1) {
                   indn -= 1;
                   p[ind] += p[indn];
               }
           }
       }
}

/*
   ! N.B. Argument NRBAR has been removed.

!   ALGORITHM AS 154.3  APPL. STATIST. (1980) VOL.29, P.311

!   FORTRAN VERSION OF REVISED VERSION OF ALGORITHM AS 75.1
!   APPL. STATIST. (1974) VOL.23, P.448
!   SEE REMARK AS R17 APPL. STATIST. (1976) VOL.25, P.323
*/


void gms::stat::inclu2(const int32_t np,
            const float weight,
            float * __restrict __attribute__((aligned(64))) xnext,
            float * __restrict __attribute__((aligned(64))) xrow,
            const float ynext,
            float * __restrict __attribute__((aligned(64))) d,
            float * __restrict __attribute__((aligned(64))) rbar,
            float * __restrict __attribute__((aligned(64))) thetab,
            float &ssqerr,
            float &recres,
            int32_t &irank,
            int32_t &ifault) {
    /*
       !   INVOKING THIS SUBROUTINE UPDATES D, RBAR, THETAB, SSQERR
!   AND IRANK BY THE INCLUSION OF XNEXT AND YNEXT WITH A
!   SPECIFIED WEIGHT. THE VALUES OF XNEXT, YNEXT AND WEIGHT WILL
!   BE CONSERVED.  THE CORRESPONDING VALUE OF RECRES IS CALCULATED.

    */
    ifault = 1;
    if(weigth<=0.0f) { return;}
    register float xk,rbthis,xk,xi,cbar,sbar;
    float wt,y,di,dpi;
    register int32_t ithisr,i1;
    y = ynext;
    wt = weight;
    if(wt<=0.0f)
#if defined __INTEL_COMPILER || defined __ICC
    __assume_aligned(xrow,64);
    __assume_aligned(xnext,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
    xrow = (float*)__builtin_assume_aligned(xrow,64);
    xnext= (float*)__builtin_assume_aligned(xnext,64);
#endif
#if defined __INTEL_COMPILER || defined __ICC
#pragma vector aligned
#pragma vector always
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(xrow:64,xnext:64)
#endif
     for(int32_t i = 0; i != np; ++i) {
         xrow[i] = xnext[i];
     } 
    recres = 0.0f;
    ifault = 0;
    ithisr = 0;
    __assume_aligned(d,64);
    __assume_aligned(rbar,64);
    __assume_aligned(thetab,64);
    for(int32_t i = 0; i != np; ++i) {
        if(xrow[i]==0.0f) {
            ithisr = ithisr + np - 1;
        }
        else {
            xi = xrow[i];
            di = d[i];
            dpi = di + wt * xi * xi;
            d[i] = dpi;
            cbar = di / dpi;
            sbar = wt * xi / dpi;
            wt = cbar * wt;
            if(i != np) {
                i1 = i + 1;
#if defined __INTEL_COMPILER || defined __ICC		
#pragma vector aligned
#pragma vector always
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(xrow:64,rbar:64)
#endif
                for(int32_t k = i1; k != np; ++i1) {
                    ithisr += 1;
                    xk = xrow[k];
                    rbthis = rbar[ithisr];
                    xrow[k] = xk - xi * rbthis;
                    rbar[ithisr] = cbar * rbthis + sbar * xk;
                }
            }
            xt = y;
            y = xk - xi * thetab[i];
            thetab[i] = cbar * thetab[i] + sbar * xk;
            if(di==0.0f) goto lab_40;
        }
    }
    ssqerr = ssqerr + wt * y * y;
    recres = y * std::sqrt(wt);
    return;
lab_40:
    irank += 1;
}

/*
   !   ALGORITHM AS 154.4  APPL. STATIST. (1980) VOL.29, P.311

!   REVISED VERSION OF ALGORITHM AS 75.4
!   APPL. STATIST. (1974) VOL.23, P.448
!   INVOKING THIS SUBROUTINE OBTAINS BETA BY BACKSUBSTITUTION
!   IN THE TRIANGULAR SYSTEM RBAR AND THETAB.
*/


void gms::stat::regres(const int32_t np,
            const int32_t nrbar,
            float * __restrict __attribute__((aligned(64))) rbar,
            float * __restrict __attribute__((aligned(64))) thetab,
            float * __restrict __attribute__((aligned(64))) beta) {
    register float bi;
    register int32_t im,jm,ithisr;
#if defined __INTEL_COMPILER || defined __ICC
    __assume_aligned(rbar,64);
    __assume_aligned(thetab,64);
    __assume_aligned(beta,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
    rbar   = (float*)__builtin_assume_aligned(rbar,64);
    thetab = (float*)__builtin_assume_aligned(thetab,64);
    beta   = (float*)__builtin_assume_aligned(beta,64);
#endif
    ithisr = nrbar;
    im = np;
    for(int32_t i = 0; i != np; ++i) {
        bi = thetab[im];
        if(im!=np) {
            i1 = i - 1;
            jm = np;
#if defined __INTEL_COMPILER
#pragma vector aligned
#pragma vector always
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(rbar:64,beta:64)
#endif
            for(int32_t j = 0; j != i1; ++j) {
                bi = bi - rbar[ithisr] * beta[jm];
                ithisr -= 1;
                jm -= 1;
            }
        }
        beta[im] = bi;
        im -= 1;
    }
}


void gms::stat::karma(const int32_t ip,
           const int32_t iq,
	   const int32_t ir,
	   float * __restrict __attribute__((aligned(64))) phi,
	   float * __restrict __attribute__((aligned(64))) theta,
	   float * __restrict __attribute__((aligned(64))) a,
	   float * __restrict __attribute__((aligned(64))) p,
	   float * __restrict __attribute__((aligned(64))) v,
	   const int32_t n,
	   float * __restrict __attribute__((aligned(64))) w,
	   float * __restrict __attribute__((aligned(64))) resid,
	   float &sumlog,
	   float &ssq,
	   int32_t iupd,
	   const float delta,
	   float * __restrict __attribute__((aligned(64))) e,
	   int32_t &nit) {

     float     wnext, a1, dt, et, ft, ut, g;
     int32_t   i, ii, ind, inde, indn, indw, ir1, j, l;
     ir1 = ir - 1;
#if defined __INTEL_COMPILER || defined __ICC
     __assume_aligned(e,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
     e = (float*)__builtin_assume_aligned(e,64);
#endif
#if defined __INTEL_COMPILER || defined __ICC
#pragma vector aligned
#pragma vector always
#elif defined __GNUC_ && !defined __INTEL_COMPILER
#pragma omp simd aligned(e:64)
#endif
     for(i = 0; i != ir; ++i) {
         e[i] = 0.0f;
     }
     inde = 1;
#if defined __INTEL_COMPILER || defined __ICC
     __assume_aligned(phi,64);
     __assume_aligned(theta,64);
     __assume_aligned(a,64);
     __assume_aligned(p,64);
     __assume_aligned(v,64);
     __assume_aligned(w,64);
     __assume_aligned(resid,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
     phi   = (float*)__builtin_assume_aligned(phi,64);
     theta = (float*)__builtin_assume_aligned(theta,64);
     a     = (float*)__builtin_assume_aligned(a,64);
     p     = (float*)__builtin_assume_aligned(p,64);
     v     = (float*)__builtin_assume_aligned(v,64);
     w     = (float*)__builtin_assume_aligned(w,64);
     resid = (float*)__builtin_assume_aligned(resid,64);
     //  FOR NON-ZERO VALUES OF NIT, PERFORM QUICK RECURSIONS.
     if(nit==0) {
        for(i = 0; i != n; ++i) {
            wnext = w[i];
	    // Prediction
	    if(iupd != 1 || i != 0) {
               dt = 0.0f;
	       if(ir != 1) dt = p[ir+1];
	       if(dt < delta) goto L100;
	       a1 = a[0];
	       if(ir != 1) {
                  for(j = 0; j != ir1; ++j) {
                      a[j] = a[j+1];
		  }
	       }
	       a[ir] = 0.0f;
	       if(ip != 0) {
#if defined __INTEL_COMPILER || defined __ICC
#pragma ivdep
#pragma vector aligned
#pragma vector always
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(a:64,phi:64)
#endif
                    for(j = 0; j != ip; ++j) {
                        a[j] = a[j] + phi[j] * a1;
		    }
	       }
	       ind = 0;
	       indn = ir;
	       for(l = 0; l != ir; ++l) {
#if defined __INTEL_COMPILER || defined __ICC
#pragma vector aligned
#pragma vector always
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(p:64,v:64)
#endif
                   for(j = l; j != ir; ++j) {
                       ind += 1;
		       p[ind] = v[ind];
		       if(j != ir) {
                          indn += 1;
			  p[ind] = p[ind] + p[indn];
		       }
		   }
	       }
	    }

	    // Updating
	    ft = p[0];
	    ut = wnext - a[0];
	    if(ir != 1) {
               ind = ir;
	       for(j = 1; j != ir; ++j) {
                   g = p[j] / ft;
		   a[j] = a[j] + g * ut;
		   for(l = j; l != ir; ++l) {
                       ind += 1;
		       p[ind] = p[ind] - g * p[l];
		   }
	       }
	    }
	    a[0] = wnext;
#if defined __INTEL_COMPILER || defined __ICC
#pragma vector aligned
#pragma vector always
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(p:64)
#endif
            for(int32_t i = 0; i != ir; ++i) {
                p[i] = 0.0f;
	    }
	    resid[i] = ut / ceph_sqrtf(rf);
	    e[inde] = resid[i];
	    inde += 1;
	    if(inde > iq) inde = 1;
	    ssq = ssq + ut * ut / ft;
	    sumlog = sumlog + ceph_logf(ft);
	}
	nit = n;
	return;
     }
     //  QUICK RECURSIONS
     i = 1;
L100:
     nit = i - 1;
     for(ii = i; ii != n; ++ii) {
         et = w[ii];
	 indw = ii;
	 if(ip != 0) {
#if defined __INTEL_COMPILER || defined __ICC
#pragma vector aligned
#pragma simd reduction(-:et)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd reduction(-:et) aligned(theta:64,e:64)
#endif
              for(j = 0; j != iq; ++j) {
                  inde -= 1;
		  if(inde==0) inde = iq;
		  et = et - theta[j] * e[inde];
	      }
	 }
	 e[inde] = et;
	 resid[ii] = et;
	 ssq = ssq + et * et;
	 inde += 1;
	 if(inde > iq) inde = 1;
     }

}





