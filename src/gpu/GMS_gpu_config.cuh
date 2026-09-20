
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

#ifndef _GMS_GPU_CONFIG_CUH_
#define _GMS_GPU_CONFIG_CUH_
	
#include <stdio.h>


#if !defined (FORTRAN_CALLABLE)
#define FORTRAN_CALLABLE 1
#else
#define FORTRAN_CALLABLE 0
#endif


#if !defined (CHECK_FORTRAN_ARRAYS)
#define CHECK_FORTRAN_ARRAYS 1
#endif

#if !defined (CHECK_FORTRAN_SCALAR_ARGS)
#define CHECK_FORTRAN_SCALAR_ARGS 1
#endif

#if defined _DEBUG
#define GMS_CUDA_DEBUG_ON 1
#else
#define GMS_CUDA_DEBUG_OFF 0
#endif

// Error handling macro
#if (GMS_CUDA_DEBUG_ON) == 1
#define GMS_CUDA_DEBUG_CHECK(func) do { \
 (status) = (func);  \
 if(cudaSuccess != (status)) { \
 fprintf(stderr, "CUDA Runtime Failure: (line %d of file %s) : \n\t" \
 "%s returned 0x%x (%s)\n", \
 __LINE__ , __FILE__ , #func,status, cudaGetErrorString(status));  \
 goto Error; \
 }    \
	} while(0) ;
#else
#define GMS_CUDA_CHECK(func) do { \
  status = (func); \
  if(cudaSuccess != (status))  { \
     goto Error;              \       
  } \
	} while(0);
#endif

#define GMS_CUBLAS_STAT_CHECK(stat)                                                                          \
    do {                                                                                           \
        cublasStatus_t err_ = (stat);                                                               \
        if (err_ != CUBLAS_STATUS_SUCCESS) {                                                       \
            printf("Abnormal cublas termination: error %d at %s:%d\n", err_, __FILE__, __LINE__);  \
            goto Error;                                                                            \
        }                                                                                          \
    } while (0)

#if !defined(SOFTWARE_PREFETCH)
#define SOFTWARE_PREFETCH 0
#endif

// Workaround for Fortran optional argument on C-side.
#if !defined (FORTRAN_OPTIONAL)
#define FORTRAN_OPTIONAL 1
#endif

#if !defined (PRINT_ERROR_TO_SCREEN)
#define PRINT_ERROR_TO_SCREEN 1
#endif

#if !defined (PRINT_ERROR_TO_FILE)
#define PRINT_ERROR_TO_FILE 1
#endif



#define GMS_CUDA_ASSERT(predicate) if ( ! (predicate)) _asm_("int3");
// helper option (not used)
#if 0
do { if(!(predicate)) {fprintf(stderr, "Asserion failed: %s at line %d in file %s\n", \
		#predicate, __LINE__,__FILE__); \
		_asm_("int 3"); \
	      } \
 } while(0); 

#endif


//
// Impotrant:
//             Set this value to '1' 
//			  if your GPU memory has more then 4 GiB.
//
#if !defined (GPU_LARGE_MEM_SPACE)
#define GPU_LARGE_MEM_SPACE 0
#endif




#if !defined (REPORT_ERROR)
#define REPORT_ERROR(msg) fprintf(stderr,"%s at line %d in file %s\n", \
	msg,__LINE__,__FILE__);
#endif

#if !defined (LOG_ACTIVITY)
#define LOG_ACTIVITY 1
#endif

#if !defined (GMS_CUDA_LOG)
#define GMS_CUDA_LOG(msg) printf("Logger: %s %s %s at line %d in file %s\n", \
	__DATE__,__TIME__ ,msg, __LINE__,__FILE__);
#endif

#if !defined (HOST_ALIGN32)
#define HOST_ALIGN32 32
#endif

#if !defined (HOST_ALIGN64) // Cache aware alignment.
#define HOST_ALIGN64 64
#endif
			       
#if !defined(PROFILE_HOST_TO_DEVICE)
#define PROFILE_HOST_TO_DEVICE 1
#endif
			       
#if !defined(PROFILE_DEVICE)
#define PROFILE_DEVICE 1
#endif


#include <stdlib.h>

//
//	Log error message to file and call exit
//

void exit_fatal(const char * msg,const char * fname) {

	assert(NULL != msg && NULL != fname);
	printf("***FATAL-ERROR***\n");
	printf(" %s\n",msg);

	FILE * fp = NULL;
	fp = fopen(fname,"a+");
	if(NULL == fp) {
		fprintf(fp, "FATAL ERROR: %s\n",msg);
		fclose(fp);
	}
	exit(EXIT_FAILURE);
}

//
// Log error message to file and call exit
// FATAL CUDA runtime error
//
void fatal_gpu_error(const char *msg, 
                     cudaError cuerr,
                     const char * fname) {
	assert(NULL != msg);
	printf("***CUDA-RUNTIME***: FATAL-ERROR\n");
	printf("%s\n",msg);
	FILE* fp = NULL;
        fp = fopen(fname,"a+");
        if(fp) {
		cudaDeviceProp dp;
		cudaError stat1,stat2;
		int dev = -1;
		stat1 = cudaGetDevice(&dev);
		stat2 = cudaGetDeviceProperties(&dp,dev);
		if(stat1 == cudaSuccess && stat2 == cudaSuccess){
		   fprintf(fp,"\tCUDA-ERROR: !!!!-- [%s] --!!!! returned by device: [%s]\n",
				                 cudaGetErrorString(cuerr),dp.name);
		   fprintf(fp, " %s \n",msg);
		}
		else {
			 fprintf(fp,"\tCUDA-ERROR: !!!-- [%s] --!!! \n",msg);
		}
		fclose(fp);
	}
	exit(EXIT_FAILURE);
}


#endif /*_GMS_GPU_CONFIG_CUH_*/
