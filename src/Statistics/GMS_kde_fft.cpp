

#include <algorithm>
#include <cmath>
#include "GMS_kde_fft.h"
#include "GMS_dyn_array.h"
#include "GMS_fft.h"

void 
gms::math 
::denest(float * __restrict__   dt,
         const int32_t          ndt,
         float                  dlo,
         float                  dhi,
         float                  window,
         float * __restrict__   ft,
         float * __restrict__   smooth,
         const int32_t          nft,
         int32_t                ical,
         int32_t * __restrict__ ifault) 
{
     
     /*
         ALGORITHM AS 176 APPL. STATIST. (1982:) VOL.314 NO.1
C
C FIND DENSITY ESTIMATE BY tKERNEL METHOD USING GAUSSIAN
C K(ERNEL. THE INTERVAL ON WHICH THE ESTIMATE IS EVALUATED
C HAS END POINTS DLO AND DHI. IF ICAL IS NOT ZERO
C THEN IT IS ASSUMED THAT THE ROUTINE HAS BEEN
C CALLED BEFORE WITH THE SAME DATA AND END POINTS
C AND THAT THE ARRAY FT HAS NOT BEEN ALTERED.


     */
    
      int32_t     ifac[15];
      float       step,ainc,hw,fac1;
      float       fnft,rj,dlo1,fndt;
      float       t0,t1,fac;
      int32_t     ii,k,nft2,j,i;
      int32_t     jj,jmax,jhi,sval; 
      int32_t     j2lo,j1,j2;
      const float zero = 0.0f;
      const float one  = 1.0f;
      const float thir2= 32.0f;
      const float big  = 30.0e+00f;
      const float atan1= 0.785398163f;
      const int32_t kftlo   = 5;
      const int32_t kfthi   = 20;
      
      if(window<zero) goto L92;
      if(dlo>=dhi)    goto L93;
      ii = 32;
      for(k=kftlo; k<=kfthi; ++k) {
          if(ii==nft) goto L2;
          ii += ii;
      }
      *ifault = 1;
      return;
L2:   
      fnft = (float)nft;
      step = (dhi-dlo)/fnft;
      fndt = (float)ndt;
      ainc = one/(fndt*step);
      nft2 = nft/2;
      hw   = window/step;
      t0   = atan1*hw/fnft;
      fac1 = thir2*t0*t0;
      if(ical!=0) goto L10;
      dlo1 = dlo-step;
      for(j=0; j<nft; ++j) ft[j] = zero;
      for(i=0; i<ndt; ++i) {
          jj = (dt[i]-dlo1)/step;
          if(jj>1 && jj<nft) ft[jj] += ainc;
      }
      //float wsave[2*nft]; // **********VLA array****************.
      gms::darray_r4_t wsave(static_cast<std::size_t>(2*nft));      
      __ogg_fdrffti(nft, &wsave.m_data[0],&ifac[0]);
      __ogg_fdrfftf(nft, &ft[0], &wsave.m_data[0],&ifac[0]);
      // call fft
L10:
      jhi       = (int32_t)std::sqrt(big/fac1);
      jmax      = std::min(nft2-1,jhi);
      smooth[0] = ft[0];
      rj        = zero;
      for(j=0; j<jmax; ++j) {
          rj         += one;
          fac        =  std::exp(-fac1*rj*rj);
          j1         = j+1;
          j2         = j1+nft2;
          smooth[j1] = fac*ft[j1];
          smooth[j2] = fac*ft[j2];
      } 
      sval = jhi+1-nft2;
      if(sval==1) {
             j2lo = jhi+2;
             for(j1=j2lo; j1<=nft2; ++j1) {
                 j2 = j1+nft2;
                 smooth[j1] = zero;
                 smooth[j2] = zero;
             }
             smooth[nft2+1] = zero;
       }
      else if(sval==2) 
             smooth[nft2+1] = zero;
      else if(sval==3)    
             smooth[nft2+1] = std::exp(-fac1*fnft*fnft)*ft[nft2+1];
       // call ifft 
       __ogg_fdrfftb(nft, &smooth[0], &wsave.m_data[0], &ifac[0]);
      for(j=0; j<nft; ++j) if(smooth[j]<zero) smooth[j] = zero;
      *ifault = 0;
       return;
L92:
      *ifault = 2;
       return;
L93:
      *ifault = 3;  
       return;
}