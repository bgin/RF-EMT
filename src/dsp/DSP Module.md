# DSP Module

This directory contains numerical DSP and mathematical kernels used by the RF-EMT project, including scalar and SIMD-optimized implementations for common trigonometric and arithmetic operations.

| File | Description |
| --- | --- |
| `GMS_32f_64f_add_64f_avx.c` | AVX-optimized implementation of a 32-bit float to 64-bit float addition kernel. |
| `GMS_32f_64f_add_64f_avx.h` | Header declarations for the AVX 32f/64f add routine. |
| `GMS_32f_64f_mul_64f_avx.c` | AVX-optimized implementation of a 32-bit float to 64-bit float multiply kernel. |
| `GMS_32f_64f_mul_64f_avx.h` | Header declarations for the AVX 32f/64f multiply routine. |
| `GMS_32f_acos_32f.c` | Scalar inverse-cosine implementation for float arrays. |
| `GMS_32f_acos_32f.h` | Public API for the scalar `acos` routine. |
| `GMS_32f_acos_32f_avx512.c` | AVX-512 optimized `acos` implementation for float vectors. |
| `GMS_32f_acos_32f_avx512.h` | Header for the AVX-512 `acos` implementation. |
| `GMS_32f_asin_32f.c` | Scalar inverse-sine implementation for float arrays. |
| `GMS_32f_asin_32f.h` | Public API for the scalar `asin` routine. |
| `GMS_32f_asin_32f_avx512.c` | AVX-512 optimized `asin` implementation for float vectors. |
| `GMS_32f_asin_32f_avx512.h` | Header for the AVX-512 `asin` implementation. |
| `GMS_32f_atan_32f.c` | Scalar inverse-tangent implementation for float arrays. |
| `GMS_32f_atan_32f.h` | Public API for the scalar `atan` routine. |
| `GMS_32f_atan_32f_avx512.c` | AVX-512 optimized `atan` implementation for float vectors. |
| `GMS_32f_atan_32f_avx512.h` | Header for the AVX-512 `atan` implementation. |
| `GMS_32f_cos_32f.c` | Scalar cosine evaluation routine for float data. |
| `GMS_32f_cos_32f.h` | Header with cosine function declarations and helpers. |
| `GMS_32f_cos_32f_avx512.c` | AVX-512 optimized cosine implementation for float vectors. |
| `GMS_32f_cos_32f_avx512.h` | Header for the AVX-512 cosine implementation. |
| `GMS_32f_exp_32f.c` | Scalar exponential function implementation for float arrays. |
| `GMS_32f_exp_32f.h` | Header for the scalar exponential routine. |
| `GMS_32f_pow_32f.c` | Scalar power-function implementation for float inputs. |
| `GMS_32f_pow_32f.h` | Header for the scalar power-function API. |
| `GMS_32f_pow_32f_avx512.c` | AVX-512 accelerated power-function implementation. |
| `GMS_32f_pow_32f_avx512.h` | Header for the AVX-512 power routine. |
| `GMS_32f_reduce_add_s32f.c` | Sum-reduction kernel for accumulating float values. |
| `GMS_32f_reduce_add_s32f.h` | Header for the reduction-add routine. |
| `GMS_32f_reduce_mul_s32f.c` | Product-reduction kernel for multiplying float values. |
| `GMS_32f_reduce_mul_s32f.h` | Header for the reduction-multiply routine. |
| `GMS_32f_sin_32f.c` | Scalar sine evaluation routine for float data. |
| `GMS_32f_sin_32f.h` | Header for the scalar sine function declarations. |
| `GMS_32f_sin_32f.hpp` | C++ inline/template implementation for sine-related functionality. |
| `GMS_32f_sin_32f_avx512.c` | AVX-512 optimized sine implementation for float vectors. |
| `GMS_32f_sin_32f_avx512.h` | Header for the AVX-512 sine implementation. |
| `GMS_32f_sin_32f_sse.c` | SSE-optimized sine implementation for float arrays. |
| `GMS_32f_sin_32f_sse.h` | Header for the SSE sine routine. |
| `GMS_32f_tan_32f.c` | Scalar tangent evaluation routine for float data. |
| `GMS_32f_tan_32f.h` | Header for the scalar tangent routine. |
| `GMS_32f_tan_32f.hpp` | C++ inline/template implementation for tangent-related functionality. |
| `GMS_IPPS_common_wrappers.h` | Common wrapper definitions for IPP-style DSP integration. |
| `GMS_IPPS_computational_wrappers.h` | Wrapper layer for computational DSP function calls and compatibility helpers. |
