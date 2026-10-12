#ifndef GMS_CABS_VEC_CUDA_H
#define GMS_CABS_VEC_CUDA_H

#include <cuda_runtime.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Host wrapper function to launch the CUDA kernels.
 * 
 * @param d_re          Device pointer to the real part array (Input)
 * @param d_im          Device pointer to the imaginary part array (Input)
 * @param d_cabs        Device pointer to the output magnitude array (Output)
 * @param n             Total number of elements in the vector
 * @param use_unrolled  Set to true to launch the kernel with #pragma unroll 4
 * @param stream        CUDA stream for asynchronous execution (pass 0 or cudaStreamDefault for default)
 */
void launch_cabsv_cuda(const float* d_re, const float* d_im, float* d_cabs, int n, bool use_unrolled = false, cudaStream_t stream = 0);

#ifdef __cplusplus
}
#endif

#endif // GMS_CABS_VEC_CUDA_H
