
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

#ifndef __GMS_ATOMICOPS_CUH__
#define __GMS_ATOMICOPS_CUH__

#define max(a,b) (a > b ? a : b)


#if defined (__CUDA_ARCH__) && __CUDA_ARCH__ < 200




__device__ void atomicFloatAdd(float* __restrict__ address, 
                               float val) 
{ 
	float old = *address, assumed; 
	do { 
		assumed = old; 
		old = __int_as_float(atomicCAS((unsigned int*)address, 
									   __float_as_int(assumed), 
									   __float_as_int(val + assumed))); 
	} while (assumed != old);
}

#endif




inline __device__ void atomicMax(float* __restrict__ address, 
                                 float val) 
{ 
	unsigned int* address_as_u = (unsigned int*) address; 
	unsigned int old = *address_as_u, assumed; 
	do
	{ 
		assumed = old; 
		old = atomicCAS(address_as_u, assumed, __float_as_int( max(val, __int_as_float(assumed) ))); 
	} while (assumed != old); 
}

inline __device__ void atomicImax(const float *   __restrict__ address, 
                                  float val, 
                                  int * __restrict__ iaddress, 
                                  int idx) 
{
	int old = *iaddress, assumed; 
	do
	{
		assumed = old;
		float themax = address[assumed];
		if (val > themax)
		{
			old = atomicCAS(iaddress, assumed, idx); 
		}
	} while (assumed != old); 
}

#endif





inline __device__ void atomicDoubleAdd(double* address, double val) 
{ 
	unsigned long long int* address_as_ull = (unsigned long long int*)address; 
	unsigned long long int old = *address_as_ull, assumed; 
	do 
	{ 
		assumed = old; 
		old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed))); 
	} while (assumed != old); 
}
inline __device__ void atomicDoubleMax(double* address, double val) 
{ 
	unsigned long long int* address_as_ull = (unsigned long long int*)address; 
	unsigned long long int old = *address_as_ull, assumed; 
	do 
	{ 
		assumed = old; 
		old = atomicCAS(address_as_ull, assumed, __double_as_longlong(max(val, __longlong_as_double(assumed)))); 
	} while (assumed != old); 
}
/*
inline __device__ double atomicDoubleAdd(double* address, double val) 
{ 
	double old = *address, assumed; 
	do { 
		assumed = old; 
		old = __longlong_as_double( atomicCAS((unsigned long long int*)address, 
											  __double_as_longlong(assumed), 
											  __double_as_longlong(val + assumed))); 
	} while (assumed != old);
}
*/
#endif /*__GMS_ATOMICOPS_CUH__*/

