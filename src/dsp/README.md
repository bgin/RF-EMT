# DSP Module

This directory contains a collection of signal-processing and numerical kernels used by the RF-EMT project. The files are organized around common vector arithmetic and math-function implementations, with several variants optimized for AVX, AVX-512, and SSE instruction sets.

| File | Description |
| --- | --- |
| `GMS_32f_64f_add_64f_avx.c` | AVX-optimized implementation of a 32-bit float to 64-bit float addition kernel. |
| `GMS_32f_64f_add_64f_avx.h` | Header declarations for the AVX 32f/64f add routines and associated constants. |
| `GMS_32f_64f_mul_64f_avx.c` | AVX-optimized implementation of a 32-bit float to 64-bit float multiply kernel. |
| `GMS_32f_64f_mul_64f_avx.h` | Header declarations for the AVX 32f/64f multiply routine. |
| `GMS_32f_acos_32f.c` | Scalar implementation of the inverse cosine function for float arrays. |
| `GMS_32f_acos_32f.h` | Public API and type definitions for the scalar `acos` implementation. |
| `GMS_32f_acos_32f_avx512.c` | AVX-512 optimized implementation of `acos` for float data. |
| `GMS_32f_acos_32f_avx512.h` | Header for the AVX-512 `acos` kernels and prototypes. |
| `GMS_32f_asin_32f.c` | Scalar implementation of the inverse sine function for float arrays. |
| `GMS_32f_asin_32f.h` | Public declarations for the scalar `asin` routine. |
| `GMS_32f_asin_32f_avx512.c` | AVX-512 optimized implementation of `asin` for float vectors. |
| `GMS_32f_asin_32f_avx512.h` | Header for the AVX-512 `asin` implementation. |
| `GMS_32f_atan_32f.c` | Scalar implementation of the inverse tangent function for float arrays. |
| `GMS_32f_atan_32f.h` | Public declarations for the scalar `atan` routine. |
| `GMS_32f_atan_32f_avx512.c` | AVX-512 optimized implementation of `atan` for float vectors. |
| `GMS_32f_atan_32f_avx512.h` | Header for the AVX-512 `atan` implementation. |
| `GMS_32f_cos_32f.c` | Scalar cosine approximation and evaluation routine for float data. |
| `GMS_32f_cos_32f.h` | Header with the cosine routine prototypes and helper declarations. |
| `GMS_32f_cos_32f_avx512.c` | AVX-512 optimized cosine kernel for high-throughput float processing. |
| `GMS_32f_cos_32f_avx512.h` | Header for the AVX-512 cosine implementation. |
| `GMS_32f_exp_32f.c` | Scalar exponential function implementation for float arrays. |
| `GMS_32f_exp_32f.h` | Public declarations for the scalar exponential routine. |
| `GMS_32f_pow_32f.c` | Scalar power-function implementation for float inputs and exponents. |
| `GMS_32f_pow_32f.h` | Header for the scalar power-function API. |
| `GMS_32f_pow_32f_avx512.c` | AVX-512 accelerated power-function implementation for float vectors. |
| `GMS_32f_pow_32f_avx512.h` | Header for the AVX-512 power routine. |
| `GMS_32f_reduce_add_s32f.c` | Sum-reduction kernel that accumulates float values into a scalar result. |
| `GMS_32f_reduce_add_s32f.h` | Header for the reduction-add routines and constants. |
| `GMS_32f_reduce_mul_s32f.c` | Product-reduction kernel for multiplying float values across a vector. |
| `GMS_32f_reduce_mul_s32f.h` | Header for the reduction-multiply routines and constants. |
| `GMS_32f_sin_32f.c` | Scalar sine approximation and evaluation routine for float data. |
| `GMS_32f_sin_32f.h` | Header for the scalar sine function declarations and support code. |
| `GMS_32f_sin_32f.hpp` | C++ inline/template implementation of sine-related functionality. |
| `GMS_32f_sin_32f_avx512.c` | AVX-512 optimized sine implementation for float vectors. |
| `GMS_32f_sin_32f_avx512.h` | Header for the AVX-512 sine kernels. |
| `GMS_32f_sin_32f_sse.c` | SSE-optimized sine implementation for float arrays. |
| `GMS_32f_sin_32f_sse.h` | Header with declarations for the SSE sine routine. |
| `GMS_32f_tan_32f.c` | Scalar tangent approximation and evaluation routine for float arrays. |
| `GMS_32f_tan_32f.h` | Header for the scalar tangent routine and related declarations. |
| `GMS_32f_tan_32f.hpp` | C++ inline/template implementation of tangent-related functionality. |
| `GMS_IPPS_common_wrappers.h` | Common wrapper definitions for IPP-style DSP function integration. |
| `GMS_IPPS_computational_wrappers.h` | High-level wrapper layer for computational DSP function calls and compatibility helpers. |

This directory focuses on numerical math and vectorized DSP primitives, with implementations tuned for different SIMD architectures and signal-processing workloads.
