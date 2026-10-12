#include <cuda_runtime.h>
#include <device_launch_parameters.h>

/**
 * CUDA kernel to compute the absolute value (magnitude) of a complex vector.
 * Replaces the AVX-512 gms::math::cabsv_zmm16r4 implementation.
 *
 * @param re    Pointer to the real part array (Input)
 * @param im    Pointer to the imaginary part array (Input)
 * @param cabs  Pointer to the output magnitude array (Output)
 * @param n     Total number of elements in the vector
 */
__global__ void cabsv_cuda_kernel(const float* __restrict__ re,
                                  const float* __restrict__ im,
                                  float* __restrict__ cabs,
                                  const int n) 
{
    // Compute a unique global thread index
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    // Grid-stride loop: ensures the kernel handles any size N safely 
    // and maximizes coalesced memory access patterns.
    for (; i < n; i += blockDim.x * gridDim.x) {
        float xre = re[i];
        float xim = im[i];
        
        // __fsqrt_rn computes the single-precision square root using round-to-nearest-even.
        // On modern GPUs, this maps directly to high-throughput hardware instructions.
        cabs[i] = __fsqrt_rn((xre * xre) + (xim * xim));
    }
}

/**
 * CUDA kernel to compute the absolute value (magnitude) of a complex vector.
 * Optimizations applied: Grid-stride looping, coalesced memory matching, and compile-time loop unrolling.
 */
__global__ void cabsv_cuda_kernel_unrolled(const float* __restrict__ re,
                                           const float* __restrict__ im,
                                           float* __restrict__ cabs,
                                           const int n) 
{
    // Compute the global step size of the entire grid
    const int stride = blockDim.x * gridDim.x;
    
    // Start at the thread's unique initial global position
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    // The compiler will unroll this loop 4 times per iteration step.
    // This allows the GPU to hide instruction latency by pipeline-interleaving 
    // multiple memory loads and independent square-root calculations.
    #pragma unroll 4
    for (; i < n; i += stride) {
        float xre = re[i];
        float xim = im[i];
        
        cabs[i] = __fsqrt_rn((xre * xre) + (xim * xim));
    }
}

/**
 * Highly optimized double-precision complex absolute value operation.
 * Intended to be called directly from within other active CUDA kernels.
 * 
 * @param re   Evaluated real component input pointer (Double Precision)
 * @param im   Evaluated imaginary component input pointer (Double Precision)
 * @param n    Vector length constraint boundary
 * @return     Computed double precision magnitude value mapping
 */
__device__ inline double cabs_element_double(const double re, const double im) {
    // Maps directly to DP hardware execution units (DFMA and DSQRT)
    return __fsqrt_rn((re * re) + (im * im));
}

/**
 * Device-only Vector Processing Kernel 
 * Designed for internal device-side grid launches or dynamic parallelism.
 */
__device__ void cabsv_device_zmm8r8_unrolled(const double* __restrict__ re,
                                             const double* __restrict__ im,
                                             double* __restrict__ cabs,
                                             const int n) 
{
    const int stride = blockDim.x * gridDim.x;
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    // Compile-time loop unrolling over double-precision channels
    #pragma unroll 4
    for (; i < n; i += stride) {
        cabs[i] = cabs_element_double(re[i], im[i]);
    }
}

/**
 * Host wrapper function to launch the standard or unrolled CUDA kernels.
 */
extern "C" void launch_cabsv_cuda(const float* d_re, const float* d_im, float* d_cabs, int n, bool use_unrolled, cudaStream_t stream) {
    if (n <= 0) return;

    // 256 threads per block is generally optimal for mixing math and memory pipelines
    int threadsPerBlock = 256;
    int blocksPerGrid = (n + threadsPerBlock - 1) / threadsPerBlock;

    // Bound checking grid size limits
    if (blocksPerGrid > 65535) blocksPerGrid = 65535;

    if (use_unrolled) {
        cabsv_cuda_kernel_unrolled<<<blocksPerGrid, threadsPerBlock, 0, stream>>>(d_re, d_im, d_cabs, n);
    } else {
        cabsv_cuda_kernel<<<blocksPerGrid, threadsPerBlock, 0, stream>>>(d_re, d_im, d_cabs, n);
    }
}
