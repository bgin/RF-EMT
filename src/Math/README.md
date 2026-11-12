# Math Module Overview

This directory contains the numerical, vector, matrix, FFT, quadrature, and SIMD math utilities used throughout RF-EMT.

| File | Description |
|---|---|
| `GMS_vector3D.h` | 3D vector type and geometry operations for position, direction, and displacement calculations. |
| `GMS_vector3D.inl` | Inline implementations for vector arithmetic, normalization, dot/cross products, and transformations. |
| `GMS_matrix2D.h` | 2D matrix abstraction and common linear-algebra helper interfaces. |
| `GMS_matrix2D.inl` | Inline 2D matrix operations, indexing, and basic numerical routines. |
| `GMS_matrix_computations.h` | High-level matrix computations such as multiplication, inverse-like operations, and decomposition support. |
| `GMS_matrix_computations.cpp` | Implementation of matrix algorithms and numerical routines used by the matrix interfaces. |
| `GMS_cartesian_system.h` | Cartesian coordinate system definitions and conversion helpers. |
| `GMS_cartesian_system.inl` | Inline cartesian transformation and coordinate manipulation utilities. |
| `GMS_coordinate_points.h` | Point and coordinate container definitions for geometric computations. |
| `GMS_coordinate_points.inl` | Inline point utilities and coordinate arithmetic. |
| `GMS_quaternion.h` | Quaternion representation and rotation math for orientation handling. |
| `GMS_quaternion.inl` | Inline quaternion arithmetic, normalization, multiplication, and interpolation routines. |
| `GMS_geometry.hpp` | Geometry kernels for spatial primitives, transforms, and physical modeling. |
| `GMS_geometry.cpp` | Core geometry implementation, including surfaces, distances, intersections, and domain operations. |
| `GMS_geocentric_coodinates.h` | Geocentric coordinate conversions and Earth-centered spatial calculations. |
| `GMS_idealized_topo.h` | Idealized terrain/topography utilities for simplified environmental models. |
| `GMS_idealized_topo.cpp` | Implementation of simplified topographic profiles and surface models. |
| `GMS_math_constants.h` | Centralized mathematical constants and physical conversion constants. |
| `GMS_fast_math_funcs.h` | Fast approximations and optimized math helper functions. |
| `GMS_fast_math_funcs.inl` | Inline implementations for optimized scalar math operations. |
| `GMS_trigfuncLUT.h` | Lookup table interface for trigonometric functions and approximation reuse. |
| `GMS_trigfuncLUT.inl` | Inline implementations of sine/cosine LUT generation and access. |
| `GMS_cephes.h` | CEPHES-based special function library for advanced transcendental functions. |
| `GMS_cephes_double.h` | Double-precision special functions derived from CEPHES-style implementations. |
| `GMS_cephes_sin_cos.h` | Sin/cos approximation and evaluation helpers. |
| `GMS_hypergeometric_func.h` | Hypergeometric function definitions and numerical evaluations. |
| `GMS_stat_distributions.h` | Statistical distributions and probability density / cumulative distribution utilities. |
| `GMS_stat_distributions.inl` | Inline distribution implementations for common analytical probability models. |
| `GMS_empirical_cdf_sampler.h` | Empirical CDF sampling interface for simulation and random data generation. |
| `GMS_gaussian_Q_approx.h` | Gaussian Q-function approximation for tail probability calculations. |
| `GMS_marcum_Q_approx.h` | Marcum Q-function approximation for radar and communication signal modeling. |
| `GMS_gaussian_Q_approx_simd.h` | SIMD-optimized Gaussian Q-function evaluation for performance-sensitive paths. |
| `GMS_marcum_Q_approx_simd.h` | SIMD-optimized Marcum Q-function routines. |
| `GMS_fft.h` | FFT API declarations and interface layer. |
| `GMS_fft.c` | FFT implementation for frequency-domain processing. |
| `GMS_real_FFT_1D.h` | Real-valued 1D FFT definitions and wrappers. |
| `GMS_real_FFT_1D.inl` | Inline code for 1D real FFT execution. |
| `GMS_real_FFT_2D.h` | Real-valued 2D FFT declarations and accessors. |
| `GMS_real_FFT_2D.inl` | Inline implementation for 2D real FFT operations. |
| `GMS_real_FFT_multiple.h` | Multi-instance real FFT utilities for batch processing. |
| `GMS_real_FFT_multiple.inl` | Inline logic for batched real FFT calculations. |
| `GMS_gauss_legendre.h` | Gauss-Legendre quadrature definitions and API surface. |
| `GMS_gauss_legendre.c` | High-accuracy Gauss-Legendre quadrature implementation. |
| `GMS_cquadpack.h` | CQUADPACK integration interface for adaptive quadrature. |
| `GMS_cquadpack.c` | Adaptive integration routines based on quadpack-style algorithms. |
| `GMS_filon_cos_quad.h` | Filon cosine quadrature API for oscillatory integrals. |
| `GMS_filon_cos_quad.cpp` | Implementation of cosine-based Filon quadrature. |
| `GMS_filon_sin_quad.h` | Filon sine quadrature API for oscillatory integrals. |
| `GMS_filon_sin_quad.cpp` | Implementation of sine-based Filon quadrature. |
| `GMS_tabulated_quadrature.h` | Interface for quadrature over tabulated or sampled data. |
| `GMS_tabulated_quadrature.cpp` | Implementation of tabulated integration routines. |
| `GMS_numerical_differentiation.h` | Numerical differentiation interfaces for approximation-based derivatives. |
| `GMS_numerical_differentiation.cpp` | Differentiation algorithms used for derivative estimation. |
| `GMS_root_finding.hpp` | Root-finding and solver library for nonlinear equations. |
| `GMS_enum_ode_integrators.hpp` | Enumeration and metadata for ODE integrator types. |
| `GMS_stochastic_rk.hpp` | Stochastic Runge-Kutta methods for SDE-like dynamics. |
| `GMS_simd_macros.h` | Cross-architecture SIMD macro definitions and portability helpers. |
| `GMS_avxvecf32.h` | AVX vectorized single-precision math utilities. |
| `GMS_avxvecf64.h` | AVX vectorized double-precision math utilities. |
| `GMS_avx512vecf32.h` | AVX-512 single-precision vector kernels. |
| `GMS_avx512vecf64.h` | AVX-512 double-precision vector kernels. |
| `GMS_avx512c8f64.h` | AVX-512 complex/packed 64-bit float operations. |
| `GMS_avxvcprocedures.h` | AVX vector conversion / procedural helper definitions. |
| `GMS_avxvcprocedures.cpp` | Implementation of AVX procedural helpers and conversion routines. |
| `GMS_avx512vcprocedures.h` | AVX-512 vector conversion and procedural support. |
| `GMS_avx512vcprocedures.cpp` | AVX-512 conversion/procedure implementation. |
| `GMS_sse_cvec4.h` | SSE-based vector container and arithmetic definitions. |
| `GMS_ref_sine_waves_LUT.h` | Reference sine-wave lookup table interfaces for waveform generation. |
| `GMS_ref_sine_waves_LUT.cpp` | Implementation of sine-wave LUT generation for signal synthesis and benchmarking. |
| `ODE_Solvers/` | Directory containing dedicated ODE/RK solver implementations and supporting routines. |
| `estimation_lib_f90/` | Fortran-based estimation library used for numerically intensive estimation routines. |

This folder is primarily a numerical math and simulation backend covering vector algebra, probability/statistics, FFTs, quadrature, differentiation, root-finding, geometry, and SIMD-accelerated kernels.
