


#include  <stdlib.h>
#include  <math.h>
#include  "GMS_nonparametric_statistics.h"
#include  "GMS_kde_fft.h"
#include  "GMS_cubint_quad.h"


void estimate_1st_moment(const float   * __restrict__ rdcyc_data, //ftab
                         const int32_t * __restrict__ xtab, // abscissas
                         float         * __restrict__ ft,  // denset (fourier-transformed data)
                         float         * __restrict__ smooth, // denset (estimated probability density)
                         const float     lo, // denset
                         const float     hi, // denset
                         const float     window,  // denset
                         const float     a, // cubint
                         const float     b, // cubint
                         const int32_t   nsamp, // number of samples must be a power of 2.
                         const int32_t   nft,   // number of FFT points
                         float * __restrict__ mean, // the 1st moment of sample.
                         int32_t  * __restrict__ cuberr, // cubint error
                         int32_t  *  __restrict denserr ) { // dnsest error
       
       float * __restrict__ integrand = NULL;
       float integral;
       float err;
       int32_t esterr;
       int32_t i;
       esterr = 9999;
       gms::math::denest(&rdcyc_data[0],nsamp,lo,hi,window,&ft[0],&smooth[0],nft,&esterr);
       *(denserr) = esterr;
       if(esterr!=0) { return;}
       integrand = (float*)malloc((size_t)nsamp*sizeof(float));
       if(integrand==NULL && nsamp!=0) { exit(EXIT_FAILURE);}
       for(i = 0; i <= nsamp; ++i) 
             integrand[i] = rdcyc_data[i]*smooth[i];
       integral = 0.0f;
       err      = nanf("");
       gms::math::cubint(ntab,&xtab[0],&integrand[0],a,b,&integral,&err);
       *(mean)   = integral;
       *(cuberr) = err;
       free(integrand);
} 
