#include <cuda_runtime.h>
#include <device_launch_parameters.h>

/**
 * ============================================================================
 * CUDA DEVICE-ONLY VECTORIZED MODIFIED BESSEL FUNCTION I0 (DOUBLE PRECISION)
 * ============================================================================
 * 
 * Ported from SSE SIMD implementation:
 * https://github.com/bgin/DSP-SIMD/blob/main/src/math_functions/GMS_bessel_i0_sse.h
 * 
 * This file implements the GPU device functions and vectorized kernels for computing 
 * the double-precision regular modified cylindrical Bessel function of order 0, I0(x).
 * 
 * It leverages the highly optimized, native mathematical hardware instructions in
 * the NVIDIA GPU streaming multiprocessors via 'cyl_bessel_i0(x)'.
 */

/**
 * Double-precision Modified Bessel Function I0 device function.
 * 
 * @param x    The input double-precision floating point value.
 * @return     The calculated value of the modified cylindrical Bessel function I0(x).
 */
__device__ __forceinline__ double gms_math_bessel_i0_r8(const double x) {
    // Maps natively to hardware-optimized polynomial evaluation on the GPU.
    // Highly accurate and performant across the entire definition domain.
    return cyl_bessel_i0(x);
}

/**
 * Device-only Vector Processing Kernel for Bessel I0.
 * Intended for device-side orchestration, dynamic parallelism, or direct 
 * invocation inside larger device pipeline routines.
 * 
 * Features:
 *  - Grid-stride looping pattern for arbitrary vector sizes.
 *  - High-performance coalesced memory reads and writes.
 *  - Compile-time loop unrolling to maximize Instruction-Level Parallelism (ILP).
 * 
 * @param input   Pointer to the input vector array on device memory (Read-only)
 * @param output  Pointer to the output vector array on device memory (Write-only)
 * @param n       Total number of elements in the vector
 */
__device__ void bessel_i0_device_zmm8r8_unrolled(const double* __restrict__ input,
                                                 double* __restrict__ output,
                                                 const int n) 
{
    // Capture the absolute global dimension stride size of the active grid execution structure
    const int stride = blockDim.x * gridDim.x;
    
    // Assign starting element relative to this specific execution thread's hardware index
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    // Instruct compiler to unroll 4 times per iteration step.
    // Pipeline-interleaves multiple memory loads and execution tracking latency.
    #pragma unroll 4
    for (; i < n; i += stride) {
        output[i] = gms_math_bessel_i0_r8(input[i]);
    }
}
