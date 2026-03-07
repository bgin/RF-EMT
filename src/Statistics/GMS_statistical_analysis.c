
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "GMS_statistical_analysis.h"
#include "GMS_nonparametric_statistics.h"
#include "GMS_descriptive_statistics.h"
#include "GMS_normality_test.h"
#include "GMS_minmax.h"





void run_stat_analysis(  const float   * __restrict__ rdcyc_data, //ftab
                         const int32_t * __restrict__ xtab, // abscissas
                         float         * __restrict__ ft,  // denset
                         float         * __restrict__ smooth, // denset
                         const float     lo, // denset
                         const float     hi, // denset
                         const float     window,  // denset
                         const float     a, // cubint
                         const float     b, // cubint
                         const int32_t   nsamp, // number of samples must be a power of 2.
                         const int32_t   nft, // number of FFT sample points;
                         const char     * __restrict__ fun_name) {
        
        float * __restrict__ sw_coef = NULL;
        float moment_1 = nanf("");
        float cuberr   = nanf("");
        float w        = 0.0f;              
        float pw       = 0.0f;                
        float srsd     = 0.0f;                  
        float svar     = 0.0f;                  
        float skew     = 0.0f;                  
        float kurt     = 0.0f;                  
        float autocor  = 0.0f;               
        float xmed     = 0.0f;     
        int32_t ifault = -1;
        int32_t esterr = -1;
        float   wlim   = 0.05f;
        int32_t n2     = nsamp/2+1   
        const bool init= false; // swilk init argument.
        sw_coef        = (float*)malloc(n2*sizeof(float));
        if(NULL==sw_coef && n2!=0) exit(EXIT_FAILURE);
        gms::math::shapiro_wilk(init,&rdcyc_data[0],nsamp,nsamp,n2,&sw_coef[0],&w,&pw,&ifault);
        if(ifault!=0) {
             printf("swilk ifault value is: %d\n",ifault);
             free(sw_coef);
             return;
        }
        printf("Normality Test [Shapiro-Wilk] results: w=%.9f,pw=%.9f\n",w,pw);
        if(pw<wlim) {
           printf("Nonparametric Statistical Analysis of the function: %s\n",fun_name); 
           estimate_1st_moment(&rdcyc_data[0],&xtab[0],&ft[0],&smooth[0],lo,hi,window,a,b,nsamp,nft,&moment_1,&cuberr,&esterr);
           printf("Estimated Mathematical Expectation (1st moment): %.9f\n",moment_1);
           printf("Input data qudrature error [cubint]: %.9f\n",cuberr);
        }
        else {
              printf("Descriptive Statistics calculations of the function: %s\n",fun_name);
	       printf("====================================================\n");
	       srsd = relsd(&rdcyc_data[0],nsamp);
	       printf("Sample Relative Standard Deviation: %.9f\n",srsd);
	       svar = var(&rdcyc_data[0],nsamp);
	       printf("Sample Variance: %.9f\n",svar);
	       skewness_kurtosis(&rdcyc_data[0],0,nsamp-1,&skew,&kurt,0);
	       printf("Skewness: %.9f, Kurtosis: %.9f\n",skew,kurt);
	       autocor = autoco(&rdcyc_data[0],nsamp);
	       printf("Autocorrelation: %.9f\n",autocor);
	       xmed    = median(&rdcyc_data[0],n2,nsamp);
	       printf("Median: %.9f\n", xmed);
	      
	}
       free(sw_coef);                
                
}               
                     
                   
         
       




