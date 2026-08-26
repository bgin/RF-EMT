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

#ifndef __GMS_CPU_PERF_TIME_SERIES_ANALYSIS_H__
#define __GMS_CPU_PERF_TIME_SERIES_ANALYSIS_H__

#include <cstdint>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <algorithm> // std::sort
#include "GMS_config.h"
#include "Timsac_iface.h"
#include "GMS_descriptive_statistics.hpp"
#include "GMS_convert_numeric_data_types.hpp"




#if !defined(DESCRIPTIVE_STATISTICS_DATA)
    #define DESCRIPTIVE_STATISTICS_DATA        \
           float  __restrict  * a    = NULL;   \
           float  __restrict  * df32 = NULL;   \
           float   w    = 0.0f;                \
           float   pw   = 0.0f;                \
           int32_t ifault = -1;                \
           float srsd = 0.0f;                  \
           float svar = 0.0f;                  \
           float skew = 0.0f;                  \
           float kurt = 0.0f;                  \
           float autocor = 0.0f;               \
           float xmid = 0.0f;                  \
           float xmean = 0.0f;                 \
           float xmidm = 0.0f;                 \
           float xmed = 0.0f;                  \
           float smin = 0.0f;                  \
           float smax = 0.0f;                  \
           float xrange = 0.0f;                \
           float xsd = 0.0f;                   \
           float xrelsd = 0.0f;                \
           float xvar = 0.0f;




/*
    Apply Time-Series analysis (Timsac) subroutine "CANARM".
    The data itself is invariant from the point of view of
    specific subroutine i.e. "CANARM".
    Attempt to calculate the descritpive statistics if result
    of Wilk-Shapiro normality test allows it. 
*/
__attribute__((hot))
__attribute__((aligned(32)))
template<int32_t len, int32_t lagh>
void
cpu_perf_time_series_canarm(const double * __restrict __attribute__((aligned(64))) data,
                                    const char   * __restrict fname,
		                    const char   * __restrict data_type,
				    const bool use_omp) {
		       
    static_assert(len <= 100000, "Input data length can not exceed -- **100000** elements!!");
     FILE * fptr = fopen(fname,"a+");
     if(__builtin_expect(NULL==fptr,0)) {
         printf("File open error: %s\n",fname);
	 std::exit(EXIT_FAILURE);
     }
    //const int32_t lagh = (int32_t)(std::sqrtf((int32_t)len));
    const int32_t len2   = len/2; // shapiro-wilk 'a' array length.
    const std::size_t  lag2len = (std::size_t)(lagh*lagh);
    const std::size_t  lag3len = (std::size_t)lag2len*len;
    constexpr float w_limit = 0.05f;
    __attribute__((aligned(64))) double acor[lagh]    = {};
    __attribute__((aligned(64))) double acov[lagh]    = {};
    __attribute__((aligned(64))) double xarcoef[lagh] = {};
    __attribute__((aligned(64))) double xv[lagh]      = {};
    __attribute__((aligned(64))) double xaic[lagh]    = {};
    __attribute__((aligned(64))) double xparcor[lagh] = {};
    __attribute__((aligned(64))) double xdicm[lagh]   = {};
    __attribute__((aligned(64))) double xb[lagh]      = {};
    __attribute__((aligned(64))) double xa[lagh]      = {};
    __attribute__((aligned(64))) int32_t xm1[lagh]    = {};
    __attribute__((aligned(64))) int32_t xm2[lagh]    = {};
    __attribute__((aligned(64))) int32_t xpo[lagh]    = {};
    double  __restrict * xw   = NULL;
    double  __restrict * xz   = NULL;
    double  __restrict * xRs  = NULL;
    double  __restrict * xchi = NULL;
    int32_t __restrict * xndt = NULL;
    double __restrict  * xdic = NULL;
    double xoaic = 0.0;
    double xmean = 0.0;
    int32_t xmo  = 0;
    int32_t xnc  = 0;
    int32_t xk   = 0;
    int32_t xl   = 0;
    DESCRIPTIVE_STATISTICS_DATA
    
    const bool init = false; // swilk init argument.
    // OpenMP multithreaded calls to _mm_malloc (using parallel sections)
    // Multithreaded allocation for large dynamic arrays.
    if(use_omp) {
    
#pragma omp parallel sections
    {
            #pragma section
	    {
               xw = reinterpret_cast<double*>(_mm_malloc(lag3len*sizeof(double),64));
	      
	    }
	    #pragma section
	    {
               xz = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	    
	    }
	    #pragma section
	    {
               xRs = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	     
	    }
	    #pragma section
	    {
              xchi = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	     
	    }
	    #pragma section
	    {
               xndt = reinterpret_cast<int32_t*>(_mm_malloc(lag2len*sizeof(int32_t),64));
	      
	    }
	    #pragma section
	    {
               xdic = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	    }
	    #pragma section
	    {
               a    = reinterpret_cast<float*>(_mm_malloc((std::size_t)len2*sizeof(float),64));
	    }
	    #pragma section
	    {
               df32 = reinterpret_cast<float*>(_mm_malloc((std::size_t)len*sizeof(float),64));
	    }
    }
      // Single thread checks the returned pointers!!
        const bool isnull =      (NULL==a)    || (NULL==xdic) || (NULL==xndt) ||
	                         (NULL==xchi) || (NULL==xRs)  || (NULL==xz)   ||
				 (NULL==xw)   || (NULL==df32);
	if(__builtin_expect(isnull,0)) {MALLOC_FAILED}
  }
  else { 

           xw = reinterpret_cast<double*>(_mm_malloc(lag3len*sizeof(double),64));
	   if(__builtin_except(NULL==xw,0)) {MALLOC_FAILED}
	   xz = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	   if(__builtin_except(NULL==xz,0)) {MALLOC_FAILED}
	   xRs = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	   if(__builtin_except(NULL==xRs,0)) {MALLOC_FAILED}
	   xchi = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	   if(__builtin_except(NULL==xchi,0)) {MALLOC_FAILED}
	   xndt = reinterpret_cast<int32_t*>(_mm_malloc(lag2len*sizeof(int32_t),64));
	   if(__builtin_except(NULL==xndt,0)) {MALLOC_FAILED}
	   xdic = reinterpret_cast<double*>(_mm_malloc(lag2len*sizeof(double),64));
	   if(__builtin_except(NULL==xdic,0)) {MALLOC_FAILED}
	   a    = reinterpret_cast<float*>(_mm_malloc((std::size_t)len2*sizeof(float),64));
	   if(__builtin_except(NULL==a,0)) {MALLOC_FAILED}
	   df32 = reinterpret_cast<float*>(_mm_malloc((std::size_t)len*sizeof(float),64));
	   if(__builtin_except(NULL==df32,0)) {MALLOC_FAILED}
    }

      
	

	autcorf_(&data[0],&len,&acov[0],&acor[0],&lagh,&xmean);
	canarmf_(&len,&lagh,&acov[0],&xarcoef[0],&lagh,&xv[0],&xaic[0],&xoaic,
	         &xmo,&xparcor[0],&xnc,&xm1[0],&xm2[0],&xw[0],&xz[0],&xRs[0],
		 &xchi[0],&xndt[0],&xdic[0],&xdicm[0],&xpo[0],&xk,&xb[0],&xl,
		 &xa[0],&lagh,&lagh);


	fprintf(fptr,"Data type: %s\n",data_type);
	fprintf(fptr,"mo=%.16f, oaic=%.16f, nc=%.16f, k=%.16f, l=%.16f\n",xmo,xoaic,xnc,xk,xl);
	fprintf(fptr, "arcoef, v, aic, parcor, dicm, b, a, m1, m2, po\n");
        for(int32_t i = 0; i != lagh; ++i) {fprintf(fptr,"%.16f %.16f %.16f %.16f %.16f %.16f %.16f %d %d %.16f\n",
         xarcoef[i],xv[i],xaic[i],xparcor[i],xdicm[i],xb[i],xa[i],xm1[i],xm2[i],xpo[i]);}
	fprintf(fptr,"w\n");
        for(int32_t i = 0; i != lag3len; ++i) {fprintf(fptr,"%.16f\n",xw[i]);}
        fprintf(fptr, "z, Rs, chi, ndt, dic\n");
        for(int32_t i = 0; i != lag2len; ++i) {fprintf(fptr, " %.16f %.16f %.16f %d %.16f\n",
           xz[i],xRs[i],xchi[i],xndt[i],xdic[i]);}
	fprintf(fptr, "End of CANARMF results dump\n");
         // Sort a samples arrays in ascending order
	//std::sort(data,data+len);
	cvrt_double_float_avx512_ptr1(&data[0],&df32[0],len);
	printf("Calling  Shapiro-Wilk normality test subroutine!!\n");
        std::sort(df32,df32+len);
        swilk(init,&df32[0],len,len,n2,&a[0],w,pw,ifault);
	if(ifault!=0) printf("swilk ifault value is: %d\n",ifault);
	fprintf(fptr,"Normality Test [Shapiro-Wilk] results: w=%.f9,pw=%.f9\n",w,pw);
	if(pw<w_limit) fprintf(fptr,"Warning!! -- 'pw' is less than normality limit -- Data is not normally distributed!!\n");
        if(pw>w_limit) {
	  
	   fprintf(fptr,"Descriptive Statistics calculations!!\n");
	   fprintf(fptr,"====================================================\n");
	   srsd = relsd(&df32[0],len);
	   fprintf(fptr,"Sample Relative Standard Deviation: %.9f\n",srsd);
	   svar = var(&df32[0],len);
	   fprintf(fptr,"Sample Variance: %.9f\n",svar);
	   skewness_kurtosis(&df32[0],0,len-1,&skew,&kurt,0);
	   fprintf(fptr,"Skewness: %.9f, Kurtosis: %.9f\n",skew,kurt);
	   autocor = autoco(&df32[0],len);
	   fprintf(fptr,"Autocorrelation: %.9f\n",autocor);
	   loc(&df32[0],len,&xmid,&xmean,&xmidm,&xmed);
	   fprintf(fptr,"Central Tendency: mid=%.9f, mean=%.9f, midm=%.9f, med=%.9f\n",
	           xmid,xmean,xmidm,xmed);
	   smin = sample_min(&df32[0],len);
	   fprintf(fptr,"Sample Min: %.9f\n",smin);
	   smax = sample_max(&df32[0],len);
	   fprintf(fptr,"Sample Max: %.9f\n",smax);
	   scale(&df32[0],len,xrange,xsd,xrelsd,xvar);
	   fprintf(fptr,"Scale Estimations: range=%.9f, sd=%.9f, relsd=%.9f, var=%.9f\n",
	           xrange,xsd,xrelsd,xvar);
	}
	fclose(fptr);
	_mm_free(df32);
	_mm_free(a);
	_mm_free(xdic);
	_mm_free(xndt);
	_mm_free(xchi);
	_mm_free(xRs);
	_mm_free(xz);
	_mm_free(xw);
}


/*
    Apply Time-Series analysis (Timsac) subroutine "MULCOR".
    The data itself is invariant from the point of view of
    specific subroutine i.e. "MULCOR".
    No descriptive statistics computations for this function.
*/
#include <string>

__attribute__((hot))
__attribute__((aligned(32)))
template<int32_t ndim, int32_t ldim, int32_t lagh>
void
cpu_perf_time_series_mulcor(const double * __restrict __attribute__((aligned(64))) mvdata, //multivariable data
                           const char   * __restrict fname,
		           const std::string * __restrict data_types){
		      
      static_assert(ndim <= 11,     "Number of dimensions can not exceed 11!!");
      static_assert(ldim <= 100000, "Number of elements per dimension can not exceed 100000!!");
      FILE * fp = fopen(fname,"a+");
      if(__builtin_expect(NULL==fp,0)) {
         printf("File open error: %s\n",fname);
	 std::exit(EXIT_FAILURE);
      }
      //const int32_t lagh = (int32_t)(2.0f*std::sqrt((float)ldim));
      const int32_t totlen = ndim*ldim;
      const std::size_t mvd_len = (std::size_t)(lagh*ndim*ndim);
      __attribute__((aligned(64))) double xmean[ndim+6];
      double * __restrict xcov = NULL;
      double * __restrict xcor = NULL;
    
      xcov = reinterpret_cast<double*>(_mm_malloc(mvd_len*sizeof(double),64));
      if(__builtin_expect(NULL==xcov,0)) {MALLOC_FAILED}
      xcor = reinterpret_cast<double*>(_mm_malloc(mvd_len*sizeof(double),64));
      if(__builtin_expect(NULL==xcor,0)) {MALLOC_FAILED}
      // Call TIMSAC MULCORF subroutine
      mulcorf_(&mvdata[0],&totlen,&ndim,&lagh,&xmean[0],&xcov[0],&xcor[0]);
    
      for(int32_t i = 0; i != ndim; ++i) {fprintf(fp,"Data types: %s\n",data_types[i].c_str());
      fprintf(fp,"Multivariate mean\n");
      for(int32_t i = 0; i != n_dim; ++i) { fprintf(fp,"%.16f\n",xmean[i]);}
      fprintf(fp1,"Multivariate Correlation and Covariance\n");
      for(int32_t i = 0; i != lagh*n_dim*n_dim; ++i) {fprintf(fp,"%.16f %.16f\n",xcor[i],xcov[i]);}
      fclose(fp1);
      _mm_free(xcor); _mm_free(xcov);
}


/*
    Apply Time-Series analysis (Timsac) subroutine "MULSPE".
    The data itself is invariant from the point of view of
    specific subroutine i.e. "MULSPE".
    No descriptive statistics computations for this function.
*/
__attribute__((hot))
__attribute__((aligned(32)))
  template<int32_t ndim, int32_t ldim, int32_t lagh>
void
cpu_perf_time_series_mulspe(const double * __restrict __attribute__((aligned(64))) mvdata, // Multidimensional data
                           const char   * __restrict fname,
		           const std::string * __restrict data_types,
			   const bool use_omp) {
       
      static_assert(ndim <= 11,     "Number of dimensions can not exceed 11!!");
      static_assert(ldim <= 100000, "Number of elements per dimension can not exceed 100000!!");
      FILE * fp = fopen(fname,"a+");
      if(__builtin_expect(NULL==fp,0)) {
         printf("File open error: %s\n",fname);
	 std::exit(EXIT_FAILURE);
      }
      //const int32_t lagh = (int32_t)(2.0f*std::sqrt((float)ldim));
      const std::size_t mvd_len = (std::size_t)(lagh*ndim*mdim);
      const int32_t totlen = ndim*ldim;
      __attribute__((aligned(64))) double xmean[ndim+6];
      __attribute__((aligned(64))) double xstat[ndim];
      // MULCOR data
      double * __restrict xcov = NULL;
      double * __restrict xcor = NULL;
      // MULSPE data
      double * __restrict xspec1 = NULL;
      double * __restrict xspec2 = NULL;
      double * __restrict xcoh1  = NULL;
      double * __restrict xcoh2  = NULL;
    
      if(use_omp) {
#pragma omp parallel sections
         {
             #pragma omp section
	        {
                   xcov =  reinterpret_cast<double*>(_mm_malloc(mvd_len*sizeof(double),64));
		}
	     #pragma omp section
	        {
                   xcor =  reinterpret_cast<double*>(_mm_malloc(mvd_len*sizeof(double),64));
		}
	     #pragma omp section
	        {
                   xspec1 =  reinterpret_cast<double*>(_mm_malloc(mvd_len*sizeof(double),64));
		}
	     #pragma omp section
	        {
                   xspec2 =  reinterpret_cast<double*>(_mm_malloc(mvd_len*sizeof(double),64));
		}
	     #pragma omp section
	        {
                   xcoh1  =  reinterpret_cast<double*>(_mm_malloc(mvd_len*sizeof(double),64));
		}
	     #pragma omp section
	        {
                   xcoh2  =  reinterpret_cast<double*>(_mm_malloc(mvd_len*sizeof(double),64)); 
	        }
	 }

	 //Single thread (main) checks for null pointers.
	 const bool isnull = (NULL==xcov)   || (NULL==xcor)  || (NULL==xspec1) ||
	                     (NULL==xspec2) || (NULL==xcoh1) || (NULL==xcoh2);
	 if(__builtin_expect(isnull,0)) {MALLOC_FAILED}
      }
      else {

            xcov =  reinterpret_cast<double*>(_mm_malloc(mvd_len*sizeof(double),64));
	    if(__bultin_expect(NULL==xcov,0)) {MALLOC_FAILED}
	    xcor =  reinterpret_cast<double*>(_mm_malloc(mvd_len*sizeof(double),64));
	    if(__builtin_expect(NULL==xcor,0)) {MALLOC_FAILED}
	    xspec1 =  reinterpret_cast<double*>(_mm_malloc(mvd_len*sizeof(double),64));
	    if(__builtin_expect(NULL==xspec1,0)) {MALLOC_FAILED}
	    xspec2 =  reinterpret_cast<double*>(_mm_malloc(mvd_len*sizeof(double),64));
	    if(__builtin_expect(NULL==xspec2,0)) {MALLOC_FAILED}
	    xcoh1  =  reinterpret_cast<double*>(_mm_malloc(mvd_len*sizeof(double),64));
	    if(__builtin_expect(NULL==xcoh1,))  {MALLOC_FAILED}
	    xcoh2  =  reinterpret_cast<double*>(_mm_malloc(mvd_len*sizeof(double),64));
	    if(__builtin_expect(NULL==xcoh2,0))  {MALLOC_FAILED}
      }
      // Call MULCORF subroutine
      mulcorf_(&mvd_data[0],&totlen,&ndim,&lagh,&xmean[0],&xcov[0],&xcor[0]);
      // Call MULSPE subroutine
      mulspef_(&tot_len,&ndim,&lagh,&lagh,&xcov[0],&xspec1[0],&xspec2[0],
               &xstat[0],&xcoh1[0],&xcoh2[0]);
     
      for(int32_t i = 0; i != ndim; ++i) {fprintf(fp,"Data types: %s\n",data_types[i].c_str());
      fprintf(fp, "Spectrum real part, imaginary part\n");
      for(int32_t i = 0; i != (int32_t)(mvd_len); ++i) { fprintf(fp,"%.16f : %.16f\n",xspec1[i],xspec2[i]);}
      fprintf(fp, "Test Statistics\n");
      for(int32_t i = 0; i != ndim; ++i) { fprintf(fp, "%.16f\n", xstat[i]);}
      fprintf(fp, "Simple coherence1, coherence2 \n");
      for(int32_t i = 0; i != (int32_t)(mvd_len); ++i) {fprintf(fp,"%.16f ,  %.16f\n",xcoh1[i],xcoh2[i]);}
      fclose(fp);
      _mm_free(xcoh2);   _mm_free(xcoh1);
      _mm_free(xspec2);  _mm_free(xspec1);
      _mm_free(xcor);    _mm_free(xcov);
}
			 


/*
    Apply Time-Series analysis (Timsac) subroutine "UNIMAR".
    The data itself is invariant from the point of view of
    specific subroutine i.e. "UNIMAR".
    Attempt to calculate the descritpive statistics if result
    of Wilk-Shapiro normality test allows it. 
*/
__attribute__((hot))
__attribute__((aligned(32)))
template<int32_t len, int32_t lagh>
void
cpu_perf_time_series_unimar(const double * __restrict __attribute__((aligned(64))) data,
                           const char   * __restrict fname,
		           const char   * __restrict data_type) {
      static_assert(len <= 1000000, "Input data can not exceed: 1000000 elements!!");
      FILE * fp = fopen(fname,"a+");
      if(__builtin_expect(NULL==fp,0)) {
         printf("File open error: %s\n",fname);
	 std::exit(EXIT_FAILURE);
      }
      const int32_t len2   = len/2; // shapiro-wilk 'a' array length.
      constexpr float w_limit = 0.05f;
      __attribute__((aligned(64))) double xv[lagh+1];
      __attribute__((aligned(64))) double xaic[lagh+1];
      __attribute__((aligned(64))) double xdaic[lagh+1];
      __attribute__((aligned(64))) double xa[lagh];
      double xmean = 0.0;
      double xvar  = 0.0;
      double xaicm = 0.0;
      double xvm   = 0.0;
      int32_t xm   = 0;
      char pad[4];
     
      DESCRIPTIVE_STATISTICS_DATA
	
      const bool init = false; // swilk init argument.
      a  = reinterpret_cast<float*>(_mm_malloc((std::size_t)len2*sizeof(float),64));
      if(__builtin_expect(NULL==a,0)) {MALLOC_FAILED}
      df32 = reinterpret_cast<float*>(_mm_malloc((std::size_t)len*sizeof(float),64));
      if(__builtin_expect(NULL==df32,0)) {MALLOC_FAILED}
      unimarf_(&data[0],&len,&lagh,&xmean,&xvar,&xv[0],&xaic[0],&xdaic[0],
               &xm,&xaicm,&xvm,&xa[0]);
     
      fprintf(fp,"Data type: %s, Method: Univariate Autoregressive AR Model Fitting\n",data_type);
      fprintf(fp,"\nmean=%.16f,var=%.16f,aicm=%.16f,vm=%.16f,xm=%d\n", xmean,
              xvar,xaicm,xvm,xm);
      fprintf(fp," V, AIC, DAIC\n");
      for(int32_t i = 0; i != lagh+1; ++i) {fprintf(fp," %.16f %.16f %.16f %.16f\n",xv[i],xaic[i],xdaic[i]);}
      fprintf(fp, "A\n");
      for(int32_t i = 0; i != lagh; ++i) {fprintf(fp," %.16f\n",xa[i]);}
      cvrt_double_float_avx512_ptr1(&data[0],&df32[0],len);
      std::sort(df32,df32+len);
      swilk(init,&df32[0],len,len,len2,&a[0],w,pw,ifault);
      if(ifault!=0) printf("swilk ifault value is: %d\n",ifault);
      fprintf(fp,"Data type: %s -- Normality Test [Shapiro-Wilk] results: w=%.f9,pw=%.f9\n",data_type,w,pw);
      if(pw<w_limit) fprintf(fp,"Warning!! -- 'pw' is less than normality limit -- Data is not normally distributed!!\n");
      if(pw>w_limit) {
	  
	   fprintf(fp,"Descriptive Statistics calculations!!\n");
	   fprintf(fp,"====================================================\n");
	   srsd = relsd(&df32[0],len);
	   fprintf(fp,"Sample Relative Standard Deviation: %.9f\n",srsd);
	   svar = var(&df32[0],len);
	   fprintf(fp,"Sample Variance: %.9f\n",svar);
	   skewness_kurtosis(&df32[0],0,len-1,&skew,&kurt,0);
	   fprintf(fp,"Skewness: %.9f, Kurtosis: %.9f\n",skew,kurt);
	   autocor = autoco(&df32[0],len);
	   fprintf(fp,"Autocorrelation: %.9f\n",autocor);
	   loc(&df32[0],len,&xmid,&xmean,&xmidm,&xmed);
	   fprintf(fp,"Central Tendency: mid=%.9f, mean=%.9f, midm=%.9f, med=%.9f\n",
	           xmid,xmean,xmidm,xmed);
	   smin = sample_min(&df32[0],len);
	   fprintf(fp,"Sample Min: %.9f\n",smin);
	   smax = sample_max(&df32[0],len);
	   fprintf(fp,"Sample Max: %.9f\n",smax);
	   scale(&df32[0],len,xrange,xsd,xrelsd,xvar);
	   fprintf(fp,"Scale Estimations: range=%.9f, sd=%.9f, relsd=%.9f, var=%.9f\n",
	           xrange,xsd,xrelsd,xvar);
	}
      fclose(fp);
      _mm_free(df32); _mm_free(a);
}
				
			



/*
    Apply Time-Series analysis (Timsac) subroutine "UNIBAR".
    The data itself is invariant from the point of view of
    specific subroutine i.e. "UNIBAR".
   
*/
__attribute__((hot))
__attribute__((aligned(32)))
template<int32_t len,int32_t lagh>
void
cpu_perf_time_series_unibar(const double * __restrict  __attribute__((aligned(64))) data,
                           const char   * __restrict fname,
		           const char   * __restrict data_type) {

      static_assert(len <= 1000000, "Input data can not exceed: 1000000 elements!!");
      FILE * fp = fopen(fname,"a+");
      if(__builtin_expect(NULL==fp,0)) {
         printf("File open error: %s\n",fname);
	 std::exit(EXIT_FAILURE);
      }
      const int32_t len2   = len/2; // shapiro-wilk 'a' array length.
      constexpr float w_limit = 0.05f;
      
      __attribute__((aligned(64))) double xv[lagh+1];
      __attribute__((aligned(64))) double xaic[lagh+1];
      __attribute__((aligned(64))) double xdaic[lagh+1];
      __attribute__((aligned(64))) double xpa[lagh];
      __attribute__((aligned(64))) double xbw[lagh+1];
      __attribute__((aligned(64))) double xsbw[lagh];
      __attribute__((aligned(64))) double xpab[lagh];
      __attribute__((aligned(64))) double xa[lagh];
      __attribute__((aligned(64))) double xpxx[128];
      double xmean = 0.0;
      double xvar  = 0.0;
      double xaicm = 0.0;
      double xvm   = 0.0;
      double xaicb = 0.0;
      double xvb   = 0.0;
      double xpn   = 0.0;
      int32_t xm   = 0;
      char pad[4];
     
      DESCRIPTIVE_STATISTICS_DATA
      
      const bool init = false; // swilk init argument.
      a  = reinterpret_cast<float*>(_mm_malloc((std::size_t)len2*sizeof(float),64));
      if(__builtin_expect(NULL==a,0)) {MALLOC_FAILED}
      df32 = reinterpret_cast<float*>(_mm_malloc((std::size_t)len*sizeof(float),64));
      if(__builtin_expect(NULL==df32,0)) {MALLOC_FAILED}
      
      unibarf_(&data[0],&len,&lagh,&xmean,&xvar,&xv[0],&xaic[0],&xdaic[0],
               &xm,&xaicm,&xvm,&xpa[0],&xbw[0],&xsbw[0],&xpab[0],&xaicb,
	       &xvb,&xpn,&xa[0],&xpxx[0]);
     
     fprintf(fp,"Data type: %s, Method: Univariate Bayesian Method of AR Model Fitting\n",data_type);
     fprintf(fp,"\nxmean=%.16f,xvar=%.16f,xaicm=%.16f,xvm=%.16f,xaicb=%.16f,xvb=%.16f,xpn=%.16f,xm=%d\n",xmean,
             xvar,xaicm,xvm,xaicb,xvb,xpn,xm);
     fprintf(fp," V, AIC, DAIC, BW\n");
     for(int32_t i = 0; i != (lagh+1); ++i) {fprintf(fp," %.16f %.16f %.16f %.16f\n",xv[i],xaic[i],xdaic[i],xbw[i]);}
     fprintf(fp, " PA, SBW, PAB, A\n");
     for(int32_t i = 0; i != lagh; ++i) {fprintf(fp," %.16f %.16f %.16f %.16f\n", xpa[i],xsbw[i],xpab[i],xa[i]);}
     fprintf(fp, " PXX\n");
     for(int32_t i = 0; i != 128; ++i) {fprintf(fp, "%.16f\n",pxx[i]);}
     cvrt_double_float_avx512_ptr1(&data[0],&df32[0],len);
     std::sort(df32,df32+len);
     swilk(init,&df32[0],len,len,len2,&a[0],w,pw,ifault);
     if(ifault!=0) printf("swilk ifault value is: %d\n",ifault);
     fprintf(fp,"Data type: %s -- Normality Test [Shapiro-Wilk] results: w=%.f9,pw=%.f9\n",data_type,w,pw);
     if(pw<w_limit) fprintf(fp,"Warning!! -- 'pw' is less than normality limit -- Data is not normally distributed!!\n");
     if(pw>w_limit) {
	  
	   fprintf(fp,"Descriptive Statistics calculations!!\n");
	   fprintf(fp,"====================================================\n");
	   srsd = relsd(&df32[0],len);
	   fprintf(fp,"Sample Relative Standard Deviation: %.9f\n",srsd);
	   svar = var(&df32[0],len);
	   fprintf(fp,"Sample Variance: %.9f\n",svar);
	   skewness_kurtosis(&df32[0],0,len-1,&skew,&kurt,0);
	   fprintf(fp,"Skewness: %.9f, Kurtosis: %.9f\n",skew,kurt);
	   autocor = autoco(&df32[0],len);
	   fprintf(fp,"Autocorrelation: %.9f\n",autocor);
	   loc(&df32[0],len,&xmid,&xmean,&xmidm,&xmed);
	   fprintf(fp,"Central Tendency: mid=%.9f, mean=%.9f, midm=%.9f, med=%.9f\n",
	           xmid,xmean,xmidm,xmed);
	   smin = sample_min(&df32[0],len);
	   fprintf(fp,"Sample Min: %.9f\n",smin);
	   smax = sample_max(&df32[0],len);
	   fprintf(fp,"Sample Max: %.9f\n",smax);
	   scale(&df32[0],len,xrange,xsd,xrelsd,xvar);
	   fprintf(fp,"Scale Estimations: range=%.9f, sd=%.9f, relsd=%.9f, var=%.9f\n",
	           xrange,xsd,xrelsd,xvar);
	}
     fclose(fp);
     _mm_free(df32); _mm_free(a);
}
				


/*
    Apply Time-Series analysis (Timsac) subroutine "EXSAR".
    The data itself is invariant from the point of view of
    specific subroutine i.e. "EXSAR".
   
*/
__attribute__((hot))
__attribute__((aligned(32)));
template<int32_t len,int32_t lagh>
void
cpu_perf_time_series_exsar( const double * __restrict  __attribute__((aligned(64))) data,
                          const char   * __restrict fname,
		          const char   * __restrict data_type) {
        static_assert(len <= 1000000, "Input data can not exceed: 1000000 elements!!");
	FILE * fp = fopen(fname,"a+");
        if(__builtin_expect(NULL==fp,0)) {
           printf("File open error: %s\n",fname);
	   std::exit(EXIT_FAILURE);
       }
	const int32_t len2   = len/2; // shapiro-wilk 'a' array length.
        constexpr float w_limit = 0.05f;
        __attribute__((aligned(64))) double xv[lagh+1];
	__attribute__((aligned(64))) double xaic[lagh+1];
	__attribute__((aligned(64))) double xdaic[lagh+1];
	__attribute__((aligned(64))) double xa1[lagh];
	__attribute__((aligned(64))) double xa2[lagh];
	double xmean = 0.0;
        double xvar  = 0.0;
        double xaicm = 0.0;
        double xsdm1 = 0.0;
        double xsdm2 = 0.0;
	char pad1[4];
        int32_t xier = 0;
        int32_t xm   = 0;
	char pad2[4];

	DESCRIPTIVE_STATISTICS_DATA

	const bool init = false; // swilk init argument.
        a  = reinterpret_cast<float*>(_mm_malloc((std::size_t)len2*sizeof(float),64));
        if(__builtin_expect(NULL==a,0)) {MALLOC_FAILED}
        df32 = reinterpret_cast<float*>(_mm_malloc((std::size_t)len*sizeof(float),64));
        if(__builtin_expect(NULL==df32,0)) {MALLOC_FAILED}
	
	exsarf_(&data[0],&len,&lagh,&xmean,&xvar,&xv[0],&xaic[0],&xdaic[0],
                &xm,&xaicm,&xsdm1,&xa1[0],&xsdm2,&xa2[0],&xier);

	 fprintf(fp,"HW Metric: %s, Maximum Likelihood Estimation\n", metric_name);
         fprintf(fp,"xmean=%.16f,xvar=%.16f,xaicm=%.16f,xsdm1=%.16f,xsdm2=%.16f,xier=%d,xm=%d\n",
                 xmean,xvar,xaicm,xsdm1,xsdm2,xier,xm);
         fprintf(fp,"V, AIC, DAIC \n");
         for(int32_t i = 0; i != lagh1; ++i) {fprintf(fp," %.16f %.16f %.16f\n", xv[i],xaic[i],xdaic[i]);}
         fprintf(fp," A1, A2 \n");
         for(int32_t i = 0; i != lagh; ++i) {fprintf(fp, " %.16f %.16f\n", xa1[i],xa2[i]);}
	 cvrt_double_float_avx512_ptr1(&data[0],&df32[0],len);
         std::sort(df32,df32+len);
         swilk(init,&df32[0],len,len,len2,&a[0],w,pw,ifault);
         if(ifault!=0) printf("swilk ifault value is: %d\n",ifault);
         fprintf(fp,"Data type: %s -- Normality Test [Shapiro-Wilk] results: w=%.f9,pw=%.f9\n",data_type,w,pw);
         if(pw<w_limit) fprintf(fp,"Warning!! -- 'pw' is less than normality limit -- Data is not normally distributed!!\n");
         if(pw>w_limit) {
	  
	   fprintf(fp,"Descriptive Statistics calculations!!\n");
	   fprintf(fp,"====================================================\n");
	   srsd = relsd(&df32[0],len);
	   fprintf(fp,"Sample Relative Standard Deviation: %.9f\n",srsd);
	   svar = var(&df32[0],len);
	   fprintf(fp,"Sample Variance: %.9f\n",svar);
	   skewness_kurtosis(&df32[0],0,len-1,&skew,&kurt,0);
	   fprintf(fp,"Skewness: %.9f, Kurtosis: %.9f\n",skew,kurt);
	   autocor = autoco(&df32[0],len);
	   fprintf(fp,"Autocorrelation: %.9f\n",autocor);
	   loc(&df32[0],len,&xmid,&xmean,&xmidm,&xmed);
	   fprintf(fp,"Central Tendency: mid=%.9f, mean=%.9f, midm=%.9f, med=%.9f\n",
	           xmid,xmean,xmidm,xmed);
	   smin = sample_min(&df32[0],len);
	   fprintf(fp,"Sample Min: %.9f\n",smin);
	   smax = sample_max(&df32[0],len);
	   fprintf(fp,"Sample Max: %.9f\n",smax);
	   scale(&df32[0],len,xrange,xsd,xrelsd,xvar);
	   fprintf(fp,"Scale Estimations: range=%.9f, sd=%.9f, relsd=%.9f, var=%.9f\n",
	           xrange,xsd,xrelsd,xvar);
	}
         fclose(fp);
	 _mm_free(df32); _mm_free(a);
}
			

/*
    Apply Time-Series analysis (Timsac) subroutine "BISPEC".
    The data itself is invariant from the point of view of
    specific subroutine i.e. "BISPEC".
  
*/
__attribute__((hot))
__attribute__((aligned(32)));
template<int32_t len,int32_t lagh>
void
cpu_perf_time_series_bispec(const double * __restrict __attribute__((aligned(64))) data,
                       const char   * __restrict fname,
		       const char   * __restrict data_type,
		       const bool use_omp) {
       static_assert(len <= 1000000, "Input data can not exceed: 1000000 elements!!");
       FILE * fp = fopen(fname,"a+");
        if(__builtin_expect(NULL==fp,0)) {
           printf("File open error: %s\n",fname);
	   std::exit(EXIT_FAILURE);
       }
       const int32_t lg12x = lagh*lagh+7;
       const std::size_t lagh_len = static_cast<std::size_t>(lg12x);\
       const int32_t len2   = len/2; // shapiro-wilk 'a' array length.
       constexpr float w_limit = 0.05f;
       const bool init = false; // swilk init argument.
       __attribute__((aligned(64))) double acor[lagh+7];
       __attribute__((aligned(64))) double acov[lagh+7];
       __attribute__((aligned(64))) double pspec1[lagh+7];
       __attribute__((aligned(64))) double psepc2[lagh+7];
       __attribute__((aligned(64))) double sig[lagh+7];
       double * __restrict mnt    = NULL;
       double * __restrict ch     = NULL;
       double * __restrict br     = NULL;
       double * __restrict bi     = NULL;
     
       double xmean = 0.0;
       double xrat = 0.0; // BISPECF result
       DESCRIPTIVE_STATISTICS_DATA
	 
       if(use_omp) {
#pragma omp parallel section
        {

       #pragma omp section
            {
              mnt = reinterpret_cast<double*>(_mm_malloc(lagh_len*sizeof(double),64));

	    }
       #pragma omp section
            {
              ch  = reinterpret_cast<double*>(_mm_malloc(lagh_len*sizeof(double),64));
	    }
       #pragma omp section
            {
              br  = reinterpret_cast<double*>(_mm_malloc(lagh_len*sizeof(double),64));
	    }
       #pragma omp section
            {
              bi  = reinterpret_cast<double*>(_mm_malloc(lagh_len*sizeof(double),64));
	    }
       #pragma omp section
	    {
               a  = reinterpret_cast<float*>(_mm_malloc((std::size_t)len2*sizeof(float),64));
	    }
       #pragma omp section
	    {
               df32 = reinterpret_cast<float*>(_mm_malloc((std::size_t)len*sizeof(float),64));
	    }

	}
	const bool isnull = (NULL==mnt) || (NULL==ch) ||
	                    (NULL==ch)  || (NULL==br) ||
	                    (NULL==bi)  || (NULL==a)  ||
			    (NULL==df32);
	if(__builtin_exppect(isnull,0)) {MALLOC_FAILED}

      }
      else {

            mnt = reinterpret_cast<double*>(_mm_malloc(lagh_len*sizeof(double),64));
	    if(__builtin_expect(NULL==mnt,0)) {MALLOC_FAILED}
	    ch  = reinterpret_cast<double*>(_mm_malloc(lagh_len*sizeof(double),64));
	    if(__builtin_expect(NULL==ch,0))  {MALLOC_FAILED}
	    br  = reinterpret_cast<double*>(_mm_malloc(lagh_len*sizeof(double),64));
	    if(__builtin_expect(NULL==br,0))  {MALLOC_FAILED}
	    bi  = reinterpret_cast<double*>(_mm_malloc(lagh_len*sizeof(double),64));
	    if(__builtin_expect(NULL==bi,0))  {MALLOC_FAILED}
	    a  = reinterpret_cast<float*>(_mm_malloc((std::size_t)len2*sizeof(float),64));
	    if(__builtin_expect(NULL==a,0))   {MALLOC_FAILED}
	    df32 = reinterpret_cast<float*>(_mm_malloc((std::size_t)len*sizeof(float),64));
	    if(__builtin_expect(NULL==df32,0)) {MALLOC_FAILED}
      }
     thirmof_(&len,&lagh,&data[0],&xmean,&acov[0],&acor[0],&mnt[0]);
     bispecf_(&len,&lagh,&data[0],&mnt[0],&pspec1[0],&pspec2[0],
              &sig[0],&br[0],&bi[0],&xrat);
    
     fprintf(fp,"Data type: %s, Bi-Spectrum Decomposition\n",data_type);
     fprintf(fp,"xrat=%.16f\n",xrat);
     fprintf(fp," %s -- Smoothed Power Spectrum-1, Power Spectrum-2 and Significance\n", metric_name);
     for(int32_t i = 0; i != lagh; ++i) { fprintf(fp, "%.16f %.16f %.16f\n", psepc1[i],pspec2[i],sig[i]);}
     fprintf(fp, " %S -- Coherence, Real part, Imaginary part\n");
     for(int32_t i = 0; i != lg12x; ++i) { fprintf(fp, "%.16f %.16f %.16f\n",ch[i],br[i],bi[i]);}
     cvrt_double_float_avx512_ptr1(&data[0],&df32[0],len);
     std::sort(df32,df32+len);
     swilk(init,&df32[0],len,len,len2,&a[0],w,pw,ifault);
     if(ifault!=0) printf("swilk ifault value is: %d\n",ifault);
     fprintf(fp,"Data type: %s -- Normality Test [Shapiro-Wilk] results: w=%.f9,pw=%.f9\n",data_type,w,pw);
     if(pw<w_limit) fprintf(fp,"Warning!! -- 'pw' is less than normality limit -- Data is not normally distributed!!\n");
     if(pw>w_limit) {
	  
	   fprintf(fp,"Descriptive Statistics calculations!!\n");
	   fprintf(fp,"====================================================\n");
	   srsd = relsd(&df32[0],len);
	   fprintf(fp,"Sample Relative Standard Deviation: %.9f\n",srsd);
	   svar = var(&df32[0],len);
	   fprintf(fp,"Sample Variance: %.9f\n",svar);
	   skewness_kurtosis(&df32[0],0,len-1,&skew,&kurt,0);
	   fprintf(fp,"Skewness: %.9f, Kurtosis: %.9f\n",skew,kurt);
	   autocor = autoco(&df32[0],len);
	   fprintf(fp,"Autocorrelation: %.9f\n",autocor);
	   loc(&df32[0],len,&xmid,&xmean,&xmidm,&xmed);
	   fprintf(fp,"Central Tendency: mid=%.9f, mean=%.9f, midm=%.9f, med=%.9f\n",
	           xmid,xmean,xmidm,xmed);
	   smin = sample_min(&df32[0],len);
	   fprintf(fp,"Sample Min: %.9f\n",smin);
	   smax = sample_max(&df32[0],len);
	   fprintf(fp,"Sample Max: %.9f\n",smax);
	   scale(&df32[0],len,xrange,xsd,xrelsd,xvar);
	   fprintf(fp,"Scale Estimations: range=%.9f, sd=%.9f, relsd=%.9f, var=%.9f\n",
	           xrange,xsd,xrelsd,xvar);
	}
     fclose(fp1);
     _mm_free(bi);   _mm_free(br);
     _mm_free(ch);   _mm_free(mnt);
     _mm_free(df32); _mm_free(a);
}


/*
    Apply Time-Series analysis (Timsac) subroutine "THIRMO".
    The data itself is invariant from the point of view of
    specific subroutine i.e. "THIRMO".
   
*/
__attribute__((hot))
__attribute__((aligned(32)));
template<int32_t len,int32_t lagh>
void
cpu_perf_time_series_thirmo(const double * __restrict __attribute__((aligned(64))) data,
                       const char   * __restrict fname,
		       const char   * __restrict data_type) {
       static_assert(len <= 1000000, "Input data can not exceed: 1000000 elements!!");
       FILE * fp = fopen(fname,"a+");
       if(__builtin_expect(NULL==fp,0)) {
         printf("File open error: %s\n",fname);
	 std::exit(EXIT_FAILURE);
       }
       const int32_t lg12x = lagh*lagh+7;
       const std::size_t lagh_len = static_cast<std::size_t>(lg12x);
       const int32_t len2   = len/2; // shapiro-wilk 'a' array length.
       constexpr float w_limit = 0.05f;
       const bool init = false; // swilk init argument.
       __attribute__((aligned(64))) double acor[lagh+7];
       __attribute__((aligned(64))) double acov[lagh+7];
       double * __restrict mnt = NULL
      
       double xmean = 0.0;
       DESCRIPTIVE_STATISTICS_DATA
       
       mnt = reinterpret_cast<double*>(_mm_malloc(lagh_len*sizeof(double),64));
       if(__builtin_expect(NULL==mnt,0)) {MALLOC_FAILED}
       a  = reinterpret_cast<float*>(_mm_malloc((std::size_t)len2*sizeof(float),64));
       if(__builtin_expect(NULL==a,0))   {MALLOC_FAILED}
       df32 = reinterpret_cast<float*>(_mm_malloc((std::size_t)len*sizeof(float),64));
       if(__builtin_expect(NULL==df32,0)) {MALLOC_FAILED}
       
       thirmof_(&len,&lagh,&data[0],&xmean,&acov[0],&acor[0],&mnt[0]);
       
        fprintf(fp,"Data type: %s Third Moments\n",data_type);
        fprintf(fp,"xmean=%.16f\n",xmean);
        fprintf(fp,"ACOV, ACOR\n");
        for(int32_t i = 0; i != lagh; ++i) { fprintf(fp, "%.16f %.16f\n", acov[i],acor[i]);}
        fprintf(fp," %S -- Third Moment\n",metric_name);
        for(int32_t i = 0; i != lg12x; ++i) { fprintf(fp, "%.16f\n",mnt[i]);}
	cvrt_double_float_avx512_ptr1(&data[0],&df32[0],len);
        std::sort(df32,df32+len);
        swilk(init,&df32[0],len,len,len2,&a[0],w,pw,ifault);
        if(ifault!=0) printf("swilk ifault value is: %d\n",ifault);
        fprintf(fp,"Data type: %s -- Normality Test [Shapiro-Wilk] results: w=%.f9,pw=%.f9\n",data_type,w,pw);
        if(pw<w_limit) fprintf(fp,"Warning!! -- 'pw' is less than normality limit -- Data is not normally distributed!!\n");
        if(pw>w_limit) {
	  
	   fprintf(fp,"Descriptive Statistics calculations!!\n");
	   fprintf(fp,"====================================================\n");
	   srsd = relsd(&df32[0],len);
	   fprintf(fp,"Sample Relative Standard Deviation: %.9f\n",srsd);
	   svar = var(&df32[0],len);
	   fprintf(fp,"Sample Variance: %.9f\n",svar);
	   skewness_kurtosis(&df32[0],0,len-1,&skew,&kurt,0);
	   fprintf(fp,"Skewness: %.9f, Kurtosis: %.9f\n",skew,kurt);
	   autocor = autoco(&df32[0],len);
	   fprintf(fp,"Autocorrelation: %.9f\n",autocor);
	   loc(&df32[0],len,&xmid,&xmean,&xmidm,&xmed);
	   fprintf(fp,"Central Tendency: mid=%.9f, mean=%.9f, midm=%.9f, med=%.9f\n",
	           xmid,xmean,xmidm,xmed);
	   smin = sample_min(&df32[0],len);
	   fprintf(fp,"Sample Min: %.9f\n",smin);
	   smax = sample_max(&df32[0],len);
	   fprintf(fp,"Sample Max: %.9f\n",smax);
	   scale(&df32[0],len,xrange,xsd,xrelsd,xvar);
	   fprintf(fp,"Scale Estimations: range=%.9f, sd=%.9f, relsd=%.9f, var=%.9f\n",
	           xrange,xsd,xrelsd,xvar);
	}
        fclose(fp);
	_mm_free(mnt);
}
				


/*
    Apply Time-Series analysis (Timsac) subroutine "AUTOCOR".
    The data itself is invariant from the point of view of
    specific subroutine i.e. "AUTOCOR".
   
*/
__attribute__((hot))
__attribute__((aligned(32)));
template<int32_t len,int32_t lagh>
void
cpu_perf_time_series_autocor(const double * __restrict __attribute__((aligned(64))) data,
                             const char   * __restrict fname,
			     const char   * __restrict data_type) {
       static_assert(len <= 1000000, "Input data can not exceed: 1000000 elements!!");
       FILE * fp = fopen(fname,"a+");
       if(__builtin_expect(NULL==fp,0)) {
         printf("File open error: %s\n",fname);
	 std::exit(EXIT_FAILURE);
       }
       const int32_t len2   = len/2; // shapiro-wilk 'a' array length.
       constexpr float w_limit = 0.05f;
       const bool init = false; // swilk init argument.
       __attribute__((aligned(64))) double acor[lagh+8];
       __attribute__((aligned(64))) double acov[lagh+8];
       double xmean = 0.0;
      
       DESCRIPTIVE_STATISTICS_DATA
       
       a  = reinterpret_cast<float*>(_mm_malloc((std::size_t)len2*sizeof(float),64));
       if(__builtin_expect(NULL==a,0))   {MALLOC_FAILED}
       df32 = reinterpret_cast<float*>(_mm_malloc((std::size_t)len*sizeof(float),64));
       if(__builtin_expect(NULL==df32,0)) {MALLOC_FAILED}
       autocorf_(&data,&len,&acov[0],&acor[0],&lagh,&xmean);
     
       fprintf(fp,"Data type: %s\n",data_type);
       fprintf(fp,"xmean=%.16f\n",xmean);
       fprintf(fp," Series Autocorrelation and Autocovariance.\n");
       for(int32_t i = 0; i != lagh; ++i) {fprintf(fp,"%.16f %.16f\n",acor[i],acov[i]);}
       cvrt_double_float_avx512_ptr1(&data[0],&df32[0],len);
       std::sort(df32,df32+len);
       swilk(init,&df32[0],len,len,len2,&a[0],w,pw,ifault);
       if(ifault!=0) printf("swilk ifault value is: %d\n",ifault);
       fprintf(fp,"Data type: %s -- Normality Test [Shapiro-Wilk] results: w=%.f9,pw=%.f9\n",data_type,w,pw);
       if(pw<w_limit) fprintf(fp,"Warning!! -- 'pw' is less than normality limit -- Data is not normally distributed!!\n");
       if(pw>w_limit) {
	  
	   fprintf(fp,"Descriptive Statistics calculations!!\n");
	   fprintf(fp,"====================================================\n");
	   srsd = relsd(&df32[0],len);
	   fprintf(fp,"Sample Relative Standard Deviation: %.9f\n",srsd);
	   svar = var(&df32[0],len);
	   fprintf(fp,"Sample Variance: %.9f\n",svar);
	   skewness_kurtosis(&df32[0],0,len-1,&skew,&kurt,0);
	   fprintf(fp,"Skewness: %.9f, Kurtosis: %.9f\n",skew,kurt);
	   autocor = autoco(&df32[0],len);
	   fprintf(fp,"Autocorrelation: %.9f\n",autocor);
	   loc(&df32[0],len,&xmid,&xmean,&xmidm,&xmed);
	   fprintf(fp,"Central Tendency: mid=%.9f, mean=%.9f, midm=%.9f, med=%.9f\n",
	           xmid,xmean,xmidm,xmed);
	   smin = sample_min(&df32[0],len);
	   fprintf(fp,"Sample Min: %.9f\n",smin);
	   smax = sample_max(&df32[0],len);
	   fprintf(fp,"Sample Max: %.9f\n",smax);
	   scale(&df32[0],len,xrange,xsd,xrelsd,xvar);
	   fprintf(fp,"Scale Estimations: range=%.9f, sd=%.9f, relsd=%.9f, var=%.9f\n",
	           xrange,xsd,xrelsd,xvar);
	}
       fclose(fp);
       _mm_free(df32); _mm_free(a);
}
				  

#endif /*__GMS_CPU_PERF_TIME_SERIES_ANALYSIS_H__*/
