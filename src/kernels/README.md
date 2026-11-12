# Kernel sources

This directory contains CPU SIMD kernels and supporting numerical routines. The `.cpp`/`.c` files generally provide implementations; matching headers declare their interfaces, while `.hpp`/`.inl` files commonly provide inline or template implementations. ISA suffixes identify variants: `sse`, `avx`/`ymm`, and `avx512`/`zmm` target progressively wider SIMD instruction sets. The `r4`/`ps` suffix denotes single precision and `r8`/`pd` denotes double precision.

| File | Contents |
|---|---|
| `BLAS-kernels/` | Subdirectory of optimized BLAS-style vector and matrix kernels. |
| `GMS_LUT_cdf_pdf_avx.cpp` | AVX implementation of lookup-table probability density and cumulative distribution functions. |
| `GMS_LUT_cdf_pdf_avx.h` | Declarations for the AVX CDF/PDF lookup-table routines. |
| `GMS_LUT_cdf_pdf_avx512.cpp` | AVX-512 implementation of lookup-table CDF/PDF routines. |
| `GMS_LUT_cdf_pdf_avx512.h` | Declarations for the AVX-512 CDF/PDF routines. |
| `GMS_LUT_cdf_pdf_sse.cpp` | SSE implementation of lookup-table CDF/PDF routines. |
| `GMS_LUT_cdf_pdf_sse.h` | Declarations for the SSE CDF/PDF routines. |
| `GMS_add_product_vec_xmm4r4.cpp` | SSE vector kernels that add products of single-precision values. |
| `GMS_add_product_vec_xmm4r4.h` | Declarations for the SSE add-product vector kernels. |
| `GMS_add_product_vec_zmm16r4.cpp` | AVX-512 single-precision add-product vector kernels. |
| `GMS_add_product_vec_zmm16r4.h` | Declarations for the AVX-512 add-product kernels. |
| `GMS_add_vec_xmm4r4.cpp` | SSE single-precision vector addition kernels. |
| `GMS_add_vec_xmm4r4.h` | Declarations for the SSE vector addition kernels. |
| `GMS_add_vec_ymm8r4.cpp` | AVX single-precision vector addition kernels. |
| `GMS_add_vec_ymm8r4.h` | Declarations for the AVX vector addition kernels. |
| `GMS_add_vec_zmm16r4.cpp` | AVX-512 single-precision vector addition kernels. |
| `GMS_add_vec_zmm16r4.h` | Declarations for the AVX-512 vector addition kernels. |
| `GMS_atan_vec_zmm16r4.cpp` | AVX-512 vectorized single-precision arctangent routines. |
| `GMS_atan_vec_zmm16r4.h` | Declarations for the AVX-512 arctangent routines. |
| `GMS_avx2_transposition_8x8.cpp` | AVX2 implementation for transposing 8×8 data blocks. |
| `GMS_avx2_transposition_8x8.h` | Declarations and support definitions for 8×8 AVX2 transposition. |
| `GMS_avx2_transposition_8x8.hpp` | Inline/template AVX2 routines for 8×8 transposition. |
| `GMS_avx512_memcpy.cpp` | AVX-512 optimized memory-copy routines. |
| `GMS_avx512_memcpy.h` | Declarations for AVX-512 memory-copy routines. |
| `GMS_avx512_memset.cpp` | AVX-512 optimized memory-fill routines. |
| `GMS_avx512_memset.h` | Declarations for AVX-512 memory-fill routines. |
| `GMS_avx512_transposition_16x16.cpp` | AVX-512 implementation for 16×16 data-block transposition. |
| `GMS_avx512_transposition_16x16.h` | Declarations for 16×16 AVX-512 transposition. |
| `GMS_avx512_transposition_4x16.cpp` | AVX-512 implementation for 4×16 data-block transposition. |
| `GMS_avx512_transposition_4x16.h` | Declarations for 4×16 AVX-512 transposition. |
| `GMS_avx512_uint64_to_double.h` | AVX-512 helpers for converting unsigned 64-bit integers to double precision. |
| `GMS_avx512_uncached_memcpy.cpp` | AVX-512 non-temporal (uncached) memory-copy routines. |
| `GMS_avx512_uncached_memcpy.h` | Declarations for AVX-512 uncached memory-copy routines. |
| `GMS_avx_deinterleave_data.h` | AVX helpers for separating interleaved data streams. |
| `GMS_avx_memcpy.cpp` | AVX optimized memory-copy routines. |
| `GMS_avx_memcpy.h` | Declarations for AVX memory-copy routines. |
| `GMS_avx_memset.cpp` | AVX optimized memory-fill routines. |
| `GMS_avx_memset.h` | Declarations for AVX memory-fill routines. |
| `GMS_avx_uncached_memcpy.cpp` | AVX non-temporal (uncached) memory-copy routines. |
| `GMS_avx_uncached_memcpy.h` | Declarations for AVX uncached memory-copy routines. |
| `GMS_bdref_avx.cpp` | AVX implementation of the BDREF radiative-transfer reference calculation. |
| `GMS_bdref_avx.h` | Declarations for the AVX BDREF routines. |
| `GMS_bdref_avx512.cpp` | AVX-512 implementation of the BDREF reference calculation. |
| `GMS_bdref_avx512.h` | Declarations for the AVX-512 BDREF routines. |
| `GMS_bdref_sse.cpp` | SSE implementation of the BDREF reference calculation. |
| `GMS_bdref_sse.h` | Declarations for the SSE BDREF routines. |
| `GMS_bessel_i0_avx.h` | AVX vector approximation of the modified Bessel function I₀. |
| `GMS_bessel_i0_avx512.h` | AVX-512 vector approximation of the modified Bessel function I₀. |
| `GMS_bessel_i0_sse.h` | SSE vector approximation of the modified Bessel function I₀. |
| `GMS_butcher_step_avx.cpp` | AVX implementation of a Butcher Runge–Kutta integration step. |
| `GMS_butcher_step_avx.h` | Declarations for the AVX Butcher integration step. |
| `GMS_butcher_step_avx512.cpp` | AVX-512 implementation of a Butcher Runge–Kutta step. |
| `GMS_butcher_step_avx512.h` | Declarations for the AVX-512 Butcher step. |
| `GMS_butcher_step_sse.cpp` | SSE implementation of a Butcher Runge–Kutta step. |
| `GMS_butcher_step_sse.h` | Declarations for the SSE Butcher step. |
| `GMS_cabs_vec_zmm16r4.cpp` | AVX-512 vectorized complex absolute-value (magnitude) operations. |
| `GMS_cabs_vec_zmm16r4.h` | Declarations for AVX-512 complex magnitude operations. |
| `GMS_cadd_vec_zmm16r4.cpp` | AVX-512 vectorized single-precision complex addition. |
| `GMS_cadd_vec_zmm16r4.h` | Declarations for AVX-512 complex addition. |
| `GMS_carg_vec_zmm16r4.cpp` | AVX-512 vectorized complex argument (phase-angle) operations. |
| `GMS_carg_vec_zmm16r4.h` | Declarations for AVX-512 complex argument operations. |
| `GMS_carithm_mean_zmm16r4.cpp` | AVX-512 arithmetic-mean calculations for complex single-precision data. |
| `GMS_carithm_mean_zmm16r4.h` | Declarations for complex arithmetic-mean routines. |
| `GMS_carithm_mean_zmm16r4.hpp` | Inline/template complex arithmetic-mean implementations. |
| `GMS_cconj_vec_zmm16r4.cpp` | AVX-512 vectorized complex conjugation. |
| `GMS_cconj_vec_zmm16r4.h` | Declarations for AVX-512 complex conjugation. |
| `GMS_ccos_vec_zmm16r4.cpp` | AVX-512 vectorized complex cosine. |
| `GMS_ccos_vec_zmm16r4.h` | Declarations for AVX-512 complex cosine. |
| `GMS_ccosh_vec_zmm16r4.cpp` | AVX-512 vectorized complex hyperbolic cosine. |
| `GMS_ccosh_vec_zmm16r4.h` | Declarations for AVX-512 complex hyperbolic cosine. |
| `GMS_cdfnorm_vec_zmm16r4.cpp` | AVX-512 vectorized normal-distribution CDF calculations. |
| `GMS_cdfnorm_vec_zmm16r4.h` | Declarations for AVX-512 normal CDF routines. |
| `GMS_cdiv_smith_vec_zmm16r4.cpp` | AVX-512 complex division using Smith's numerically stable method. |
| `GMS_cdiv_smith_vec_zmm16r4.h` | Declarations for Smith-method AVX-512 complex division. |
| `GMS_cdiv_vec_zmm16r4.cpp` | AVX-512 vectorized complex division. |
| `GMS_cdiv_vec_zmm16r4.h` | Declarations for AVX-512 complex division. |
| `GMS_ceq_vec_zmm16r4.cpp` | AVX-512 vectorized complex equality comparisons. |
| `GMS_ceq_vec_zmm16r4.h` | Declarations for complex equality comparisons. |
| `GMS_cexp_vec_zmm16r4.cpp` | AVX-512 vectorized complex exponential. |
| `GMS_cexp_vec_zmm16r4.h` | Declarations for AVX-512 complex exponential. |
| `GMS_cgemv_t.c` | Complex general matrix-vector multiplication with transposed matrix access. |
| `GMS_cgemv_t.h` | Declarations for the complex transposed matrix-vector kernel. |
| `GMS_cgt_vec_zmm16r4.cpp` | AVX-512 vectorized complex greater-than comparisons. |
| `GMS_cgt_vec_zmm16r4.h` | Declarations for complex greater-than comparisons. |
| `GMS_cgt_vec_zmm16r4.hpp` | Inline/template complex greater-than comparison implementations. |
| `GMS_clog_vec_zmm16r4.cpp` | AVX-512 vectorized complex natural logarithm. |
| `GMS_clog_vec_zmm16r4.h` | Declarations for AVX-512 complex logarithm. |
| `GMS_clt_vec_zmm16r4.cpp` | AVX-512 vectorized complex less-than comparisons. |
| `GMS_clt_vec_zmm16r4.h` | Declarations for complex less-than comparisons. |
| `GMS_cmean_cprod_zmm16r4.cpp` | AVX-512 complex mean and product-related kernels. |
| `GMS_cmean_cprod_zmm16r4.h` | Declarations for complex mean/product routines. |
| `GMS_cmean_quot_zmm16r4.cpp` | AVX-512 complex mean and quotient-related kernels. |
| `GMS_cmean_quot_zmm16r4.h` | Declarations for complex mean/quotient routines. |
| `GMS_cmul_vec_zmm16r4.cpp` | AVX-512 vectorized single-precision complex multiplication. |
| `GMS_cmul_vec_zmm16r4.h` | Declarations for AVX-512 complex multiplication. |
| `GMS_cmul_vec_zmm16r4.hpp` | Inline/template complex multiplication implementations. |
| `GMS_cneq_vec_zmm16r4.cpp` | AVX-512 vectorized complex inequality comparisons. |
| `GMS_cneq_vec_zmm16r4.h` | Declarations for complex inequality comparisons. |
| `GMS_cnorm_prod_zmm16r4.cpp` | AVX-512 complex norm and product computations. |
| `GMS_cnorm_prod_zmm16r4.h` | Declarations for complex norm/product routines. |
| `GMS_complex_common_ymm4r8.h` | Shared complex-number helpers for AVX double-precision vectors. |
| `GMS_complex_common_ymm8r4.h` | Shared complex-number helpers for AVX single-precision vectors. |
| `GMS_complex_common_zmm16r4.h` | Shared complex-number helpers for AVX-512 single-precision vectors. |
| `GMS_complex_common_zmm8r8.h` | Shared complex-number helpers for AVX-512 double-precision vectors. |
| `GMS_complex_vec1D_ymm4r8.cpp` | AVX double-precision complex 1D vector operations. |
| `GMS_complex_vec1D_ymm4r8.h` | Declarations for AVX double-precision complex-vector operations. |
| `GMS_complex_vec1D_ymm8r4.cpp` | AVX single-precision complex 1D vector operations. |
| `GMS_complex_vec1D_ymm8r4.h` | Declarations for AVX single-precision complex-vector operations. |
| `GMS_complex_vec1D_zmm16r4.cpp` | AVX-512 single-precision complex 1D vector operations. |
| `GMS_complex_vec1D_zmm16r4.h` | Declarations for AVX-512 single-precision complex-vector operations. |
| `GMS_complex_vec1D_zmm8r8.cpp` | AVX-512 double-precision complex 1D vector operations. |
| `GMS_complex_vec1D_zmm8r8.h` | Declarations for AVX-512 double-precision complex-vector operations. |
| `GMS_complex_xmm2r8.h` | SSE complex-number operations for double-precision lanes. |
| `GMS_complex_xmm4r4.h` | SSE complex-number operations for single-precision lanes. |
| `GMS_complex_ymm4r8.h` | AVX complex-number operations for double precision. |
| `GMS_complex_ymm8r4.h` | AVX complex-number operations for single precision. |
| `GMS_complex_zmm16r4.hpp` | AVX-512 single-precision complex-number operations. |
| `GMS_complex_zmm8r8.hpp` | AVX-512 double-precision complex-number operations. |
| `GMS_convert_int_to_float_avx512.c` | AVX-512 integer-to-single-precision conversion routines. |
| `GMS_convert_int_to_float_avx512.h` | Declarations for AVX-512 integer-to-float conversions. |
| `GMS_coord_conv_avx.cpp` | AVX coordinate-system conversion kernels. |
| `GMS_coord_conv_avx.h` | Declarations for AVX coordinate conversions. |
| `GMS_coord_conv_avx512.cpp` | AVX-512 coordinate-system conversion kernels. |
| `GMS_coord_conv_avx512.h` | Declarations for AVX-512 coordinate conversions. |
| `GMS_coord_conv_sse.cpp` | SSE coordinate-system conversion kernels. |
| `GMS_coord_conv_sse.h` | Declarations for SSE coordinate conversions. |
| `GMS_cos_vec_zmm16r4.cpp` | AVX-512 vectorized single-precision cosine. |
| `GMS_cos_vec_zmm16r4.h` | Declarations for AVX-512 cosine routines. |
| `GMS_cosh_vec_zmm16r4.cpp` | AVX-512 vectorized single-precision hyperbolic cosine. |
| `GMS_cosh_vec_zmm16r4.h` | Declarations for AVX-512 hyperbolic cosine routines. |
| `GMS_cpolar_zmm16r4.cpp` | AVX-512 conversion between complex rectangular and polar forms. |
| `GMS_cpolar_zmm16r4.h` | Declarations for complex polar-coordinate conversion. |
| `GMS_cpow_vec_zmm16r4.cpp` | AVX-512 vectorized complex exponentiation/power operations. |
| `GMS_cpow_vec_zmm16r4.h` | Declarations for AVX-512 complex power routines. |
| `GMS_cscal.cpp` | Complex-vector scaling kernel. |
| `GMS_cscal.h` | Declarations for complex-vector scaling. |
| `GMS_csin_vec_zmm16r4.cpp` | AVX-512 vectorized complex sine. |
| `GMS_csin_vec_zmm16r4.h` | Declarations for AVX-512 complex sine. |
| `GMS_csinh_vec_zmm16r4.cpp` | AVX-512 vectorized complex hyperbolic sine. |
| `GMS_csinh_vec_zmm16r4.h` | Declarations for AVX-512 complex hyperbolic sine. |
| `GMS_csinh_vec_zmm16r4.hpp` | Inline/template complex hyperbolic-sine implementations. |
| `GMS_csinhv_svml_zmm16r4.cpp` | AVX-512 complex hyperbolic-sine implementation using Intel SVML. |
| `GMS_csinhv_svml_zmm16r4.h` | Declarations for the SVML-based complex sinh routines. |
| `GMS_csqrt_vec_zmm16r4.cpp` | AVX-512 vectorized complex square root. |
| `GMS_csqrt_vec_zmm16r4.h` | Declarations for AVX-512 complex square root. |
| `GMS_csub_vec_zmm16r4.cpp` | AVX-512 vectorized single-precision complex subtraction. |
| `GMS_csub_vec_zmm16r4.h` | Declarations for AVX-512 complex subtraction. |
| `GMS_ctan_vec_zmm16r4.cpp` | AVX-512 vectorized complex tangent. |
| `GMS_ctan_vec_zmm16r4.h` | Declarations for AVX-512 complex tangent. |
| `GMS_ctanh_vec_zmm16r4.cpp` | AVX-512 vectorized complex hyperbolic tangent. |
| `GMS_ctanh_vec_zmm16r4.h` | Declarations for AVX-512 complex hyperbolic tangent. |
| `GMS_dd.h` | Double-double arithmetic types and helper operations. |
| `GMS_derivative_5P4P_avx.hpp` | AVX five-point/fourth-order finite-difference derivative templates. |
| `GMS_derivative_avx.cpp` | AVX numerical derivative kernels. |
| `GMS_derivative_avx.h` | Declarations for AVX derivative kernels. |
| `GMS_derivative_avx512.cpp` | AVX-512 numerical derivative kernels. |
| `GMS_derivative_avx512.h` | Declarations for AVX-512 derivative kernels. |
| `GMS_derivative_sse.cpp` | SSE numerical derivative kernels. |
| `GMS_derivative_sse.h` | Declarations for SSE derivative kernels. |
| `GMS_df.h` | Distribution-function support types and helpers. |
| `GMS_dryden_psd_avx.cpp` | AVX Dryden turbulence power-spectral-density calculations. |
| `GMS_dryden_psd_avx.h` | Declarations for AVX Dryden PSD routines. |
| `GMS_dryden_psd_avx512.cpp` | AVX-512 Dryden turbulence PSD calculations. |
| `GMS_dryden_psd_avx512.h` | Declarations for AVX-512 Dryden PSD routines. |
| `GMS_dryden_psd_sse.cpp` | SSE Dryden turbulence PSD calculations. |
| `GMS_dryden_psd_sse.h` | Declarations for SSE Dryden PSD routines. |
| `GMS_erf_vec_zmm16r4.cpp` | AVX-512 vectorized error function (erf). |
| `GMS_erf_vec_zmm16r4.h` | Declarations for AVX-512 erf routines. |
| `GMS_erfc_vec_zmm16r4.cpp` | AVX-512 vectorized complementary error function (erfc). |
| `GMS_erfc_vec_zmm16r4.h` | Declarations for AVX-512 erfc routines. |
| `GMS_erfcinv_vec_zmm16r4.cpp` | AVX-512 vectorized inverse complementary error function. |
| `GMS_erfcinv_vec_zmm16r4.h` | Declarations for AVX-512 inverse-erfc routines. |
| `GMS_euler_step_avx.hpp` | Inline AVX Euler-method ODE integration step. |
| `GMS_euler_step_avx512.hpp` | Inline AVX-512 Euler-method ODE integration step. |
| `GMS_exp_vec_zmm16r4.cpp` | AVX-512 vectorized single-precision exponential. |
| `GMS_exp_vec_zmm16r4.h` | Declarations for AVX-512 exponential routines. |
| `GMS_fast_simd_funcs_approx.h` | Fast approximate SIMD math-function helpers. |
| `GMS_fdiff1246_avx512.hpp` | AVX-512 finite-difference derivative formulas of orders 1, 2, 4, and 6. |
| `GMS_fdiff_avx512.hpp` | AVX-512 finite-difference derivative templates. |
| `GMS_geodesy_ymm4r8.cpp` | AVX double-precision geodesy computations. |
| `GMS_geodesy_ymm4r8.h` | Declarations for AVX double-precision geodesy routines. |
| `GMS_geodesy_ymm8r4.cpp` | AVX single-precision geodesy computations. |
| `GMS_geodesy_ymm8r4.h` | Declarations for AVX single-precision geodesy routines. |
| `GMS_geodesy_zmm16r4.cpp` | AVX-512 single-precision geodesy computations. |
| `GMS_geodesy_zmm16r4.h` | Declarations for AVX-512 single-precision geodesy routines. |
| `GMS_geodesy_zmm8r8.cpp` | AVX-512 double-precision geodesy computations. |
| `GMS_geodesy_zmm8r8.h` | Declarations for AVX-512 double-precision geodesy routines. |
| `GMS_helperavx2.h` | Shared AVX2 intrinsics and low-level helper routines. |
| `GMS_helperavx512f.h` | Shared AVX-512 Foundation intrinsics and helper routines. |
| `GMS_heuns_step_avx.cpp` | AVX implementation of Heun's ODE integration step. |
| `GMS_heuns_step_avx.h` | Declarations for the AVX Heun step. |
| `GMS_heuns_step_avx512.cpp` | AVX-512 implementation of Heun's ODE integration step. |
| `GMS_heuns_step_avx512.h` | Declarations for the AVX-512 Heun step. |
| `GMS_heuns_step_sse.cpp` | SSE implementation of Heun's ODE integration step. |
| `GMS_heuns_step_sse.h` | Declarations for the SSE Heun step. |
| `GMS_mat8x8_det_avx.cpp` | AVX determinant calculation for 8×8 matrices. |
| `GMS_mat8x8_det_avx.h` | Declarations for the AVX 8×8 matrix determinant routine. |
| `GMS_mueller_calculus_avx512.hpp` | AVX-512 Mueller-matrix and polarization-calculus operations. |
| `GMS_mueller_types_avx512.hpp` | AVX-512 types and data structures for Mueller calculus. |
| `GMS_nystrom_step_avx.cpp` | AVX implementation of a Nyström ODE integration step. |
| `GMS_nystrom_step_avx.h` | Declarations for the AVX Nyström step. |
| `GMS_nystrom_step_avx512.cpp` | AVX-512 implementation of a Nyström ODE integration step. |
| `GMS_nystrom_step_avx512.h` | Declarations for the AVX-512 Nyström step. |
| `GMS_nystrom_step_sse.cpp` | SSE implementation of a Nyström ODE integration step. |
| `GMS_nystrom_step_sse.h` | Declarations for the SSE Nyström step. |
| `GMS_pdf_cdf_avx.hpp` | AVX probability-density and cumulative-distribution function kernels. |
| `GMS_pdf_cdf_avx512.hpp` | AVX-512 probability-density and cumulative-distribution function kernels. |
| `GMS_pdf_cdf_sse.hpp` | SSE probability-density and cumulative-distribution function kernels. |
| `GMS_pdf_kernels_sse.h` | SSE probability-density function kernel declarations and implementations. |
| `GMS_phase_pdf_avx.cpp` | AVX phase probability-density calculations. |
| `GMS_phase_pdf_avx.h` | Declarations for AVX phase PDF routines. |
| `GMS_phase_pdf_avx512.cpp` | AVX-512 phase probability-density calculations. |
| `GMS_phase_pdf_avx512.h` | Declarations for AVX-512 phase PDF routines. |
| `GMS_phase_pdf_sse.cpp` | SSE phase probability-density calculations. |
| `GMS_phase_pdf_sse.h` | Declarations for SSE phase PDF routines. |
| `GMS_pos_to_state_avx.hpp` | AVX position-to-state conversion routines. |
| `GMS_pos_to_state_avx512.hpp` | AVX-512 position-to-state conversion routines. |
| `GMS_pos_to_state_sse.hpp` | SSE position-to-state conversion routines. |
| `GMS_q4x16_to_rmat9x16_looped_avx512.cpp` | AVX-512 conversion of batches of 4D quaternions to 3×3 rotation matrices. |
| `GMS_q4x16_to_rmat9x16_looped_avx512.h` | Declarations for batched quaternion-to-rotation-matrix conversion. |
| `GMS_rand_rmat9x16_looped_avx512.cpp` | AVX-512 generation of batches of random 3×3 rotation matrices (16-wide). |
| `GMS_rand_rmat9x16_looped_avx512.h` | Declarations for 16-wide random rotation-matrix generation. |
| `GMS_rand_rmat9x8_looped_avx512.cpp` | AVX-512 generation of batches of random 3×3 rotation matrices (8-wide). |
| `GMS_rand_rmat9x8_looped_avx512.h` | Declarations for 8-wide random rotation-matrix generation. |
| `GMS_range_rate_avx.cpp` | AVX radar range and range-rate calculations. |
| `GMS_range_rate_avx.h` | Declarations for AVX range/range-rate routines. |
| `GMS_range_rate_avx512.cpp` | AVX-512 radar range and range-rate calculations. |
| `GMS_range_rate_avx512.h` | Declarations for AVX-512 range/range-rate routines. |
| `GMS_range_rate_sse.cpp` | SSE radar range and range-rate calculations. |
| `GMS_range_rate_sse.h` | Declarations for SSE range/range-rate routines. |
| `GMS_rk38_step_avx.cpp` | AVX implementation of the Runge–Kutta 3/8 integration step. |
| `GMS_rk38_step_avx.h` | Declarations for the AVX RK 3/8 step. |
| `GMS_rk38_step_avx.hpp` | Inline/template AVX RK 3/8 integration implementation. |
| `GMS_rk3_step_avx.cpp` | AVX third-order Runge–Kutta integration step. |
| `GMS_rk3_step_avx.h` | Declarations for the AVX third-order RK step. |
| `GMS_rk3_step_avx.hpp` | Inline/template AVX third-order RK implementation. |
| `GMS_rk3v2_step_avx.cpp` | AVX implementation of a second variant of the third-order RK step. |
| `GMS_rk3v2_step_avx.h` | Declarations for the AVX RK3 variant-2 step. |
| `GMS_rk3v2_step_avx512.cpp` | AVX-512 implementation of the RK3 variant-2 step. |
| `GMS_rk3v2_step_avx512.h` | Declarations for the AVX-512 RK3 variant-2 step. |
| `GMS_rk3v2_step_avx512.hpp` | Inline/template AVX-512 RK3 variant-2 implementation. |
| `GMS_rk4_step_avx.cpp` | AVX fourth-order Runge–Kutta integration step. |
| `GMS_rk4_step_avx.h` | Declarations for the AVX fourth-order RK step. |
| `GMS_rk4_step_avx512.cpp` | AVX-512 fourth-order Runge–Kutta integration step. |
| `GMS_rk4_step_avx512.h` | Declarations for the AVX-512 fourth-order RK step. |
| `GMS_rk4_step_avx512.hpp` | Inline/template AVX-512 fourth-order RK implementation. |
| `GMS_rkg_step_avx.hpp` | Inline AVX implementation of a Runge–Kutta–Gill integration step. |
| `GMS_rkg_step_avx512.cpp` | AVX-512 Runge–Kutta–Gill integration step. |
| `GMS_rkg_step_avx512.h` | Declarations for the AVX-512 RK–Gill step. |
| `GMS_rkr2_step_avx.hpp` | Inline AVX second-order Runge–Kutta integration step. |
| `GMS_rkr2_step_avx512.hpp` | Inline AVX-512 second-order Runge–Kutta integration step. |
| `GMS_rkr4_step_avx.hpp` | Inline AVX fourth-order Runge–Kutta integration step. |
| `GMS_rkr4_step_sse.cpp` | SSE fourth-order Runge–Kutta integration step. |
| `GMS_rkr4_step_sse.h` | Declarations for the SSE fourth-order RK step. |
| `GMS_rotation_kernels_avx512.cpp` | AVX-512 vector and matrix rotation kernels. |
| `GMS_rotation_kernels_avx512.h` | Declarations for AVX-512 rotation kernels. |
| `GMS_rotation_kernels_avx512.hpp` | Inline/template AVX-512 rotation implementations. |
| `GMS_rotations_avx2_helpers.h` | AVX2 helper operations for rotations. |
| `GMS_rotations_avx512_helpers.hpp` | AVX-512 helper operations for rotations. |
| `GMS_rotations_sse_helpers.h` | SSE helper operations for rotations. |
| `GMS_setv_avx512_unroll16x.cpp` | AVX-512 vector initialization/fill routines with 16× loop unrolling. |
| `GMS_setv_avx512_unroll16x.h` | Declarations for AVX-512 16×-unrolled vector set routines. |
| `GMS_simd_cephes_cos.h` | SIMD Cephes-based cosine approximation routines. |
| `GMS_simd_cephes_sin.h` | SIMD Cephes-based sine approximation routines. |
| `GMS_simd_macros.h` | Common macros for SIMD kernel implementation. |
| `GMS_simd_utils.h` | General SIMD utility functions and helpers. |
| `GMS_simd_utils_avx512_float.h` | AVX-512 single-precision utility routines. |
| `GMS_simd_utils_constants.h` | Constants used by SIMD numerical kernels. |
| `GMS_sin_vec_zmm16r4.cpp` | AVX-512 vectorized single-precision sine. |
| `GMS_sin_vec_zmm16r4.h` | Declarations for AVX-512 sine routines. |
| `GMS_sinh_vec_zmm16r4.cpp` | AVX-512 vectorized single-precision hyperbolic sine. |
| `GMS_sinh_vec_zmm16r4.h` | Declarations for AVX-512 hyperbolic-sine routines. |
| `GMS_sleefsimddp.hpp` | SIMD double-precision elementary math routines based on SLEEF. |
| `GMS_sleefsimdsp.cpp` | SIMD single-precision elementary math routines based on SLEEF. |
| `GMS_sleefsimdsp.h` | Declarations for the SLEEF single-precision SIMD routines. |
| `GMS_spec_funcs_LUT.h` | Lookup tables for special mathematical functions. |
| `GMS_spec_funcs_zmm8r8.cpp` | AVX-512 double-precision special-function kernels. |
| `GMS_spec_funcs_zmm8r8.h` | Declarations for AVX-512 double-precision special functions. |
| `GMS_spher_grad_avx.cpp` | AVX spherical-coordinate gradient calculations. |
| `GMS_spher_grad_avx.h` | Declarations for AVX spherical-gradient routines. |
| `GMS_spher_grad_avx512.cpp` | AVX-512 spherical-coordinate gradient calculations. |
| `GMS_spher_grad_avx512.h` | Declarations for AVX-512 spherical-gradient routines. |
| `GMS_spher_grad_sse.cpp` | SSE spherical-coordinate gradient calculations. |
| `GMS_spher_grad_sse.h` | Declarations for SSE spherical-gradient routines. |
| `GMS_sse_matrix4x4.h` | SSE 4×4 matrix types and operation declarations. |
| `GMS_sse_matrix4x4.inl` | Inline SSE 4×4 matrix operations. |
| `GMS_sse_memcpy.cpp` | SSE optimized memory-copy routines. |
| `GMS_sse_memcpy.h` | Declarations for SSE memory-copy routines. |
| `GMS_sse_memset.cpp` | SSE optimized memory-fill routines. |
| `GMS_sse_memset.h` | Declarations for SSE memory-fill routines. |
| `GMS_sse_uncached_memcpy.cpp` | SSE non-temporal (uncached) memory-copy routines. |
| `GMS_sse_uncached_memcpy.h` | Declarations for SSE uncached memory-copy routines. |
| `GMS_static_cvec1D_ymm4r8.h` | Static 1D complex-vector container for AVX double precision. |
| `GMS_static_cvec1D_ymm8r4.h` | Static 1D complex-vector container for AVX single precision. |
| `GMS_static_cvec1D_zmm16r4.h` | Static 1D complex-vector container for AVX-512 single precision. |
| `GMS_static_cvec1D_zmm8r8.h` | Static 1D complex-vector container for AVX-512 double precision. |
| `GMS_stochastic_rk_avx2_pd.cpp` | AVX2 double-precision stochastic Runge–Kutta routines. |
| `GMS_stochastic_rk_avx2_pd.h` | Declarations for AVX2 double-precision stochastic RK routines. |
| `GMS_stochastic_rk_avx512.hpp` | Inline/template AVX-512 stochastic Runge–Kutta routines. |
| `GMS_stochastic_rk_avx512_pd.cpp` | AVX-512 double-precision stochastic Runge–Kutta implementation. |
| `GMS_stochastic_rk_avx512_pd.h` | Declarations for double-precision AVX-512 stochastic RK. |
| `GMS_stochastic_rk_avx512_ps.cpp` | AVX-512 single-precision stochastic Runge–Kutta implementation. |
| `GMS_stochastic_rk_avx512_ps.h` | Declarations for single-precision AVX-512 stochastic RK. |
| `GMS_swap_avx_unroll8x.cpp` | AVX vector-swap routines with 8× loop unrolling. |
| `GMS_swap_avx_unroll8x.h` | Declarations for AVX 8×-unrolled swap routines. |
| `GMS_vec_cephes_cos_ps.h` | Cephes-based SIMD cosine approximation for packed single-precision values. |
| `GMS_vecspecfuncs.h` | SIMD implementations of vectorized special mathematical functions. |
| `GMS_verner_step_avx.cpp` | AVX Verner Runge–Kutta integration step. |
| `GMS_verner_step_avx.h` | Declarations for the AVX Verner integration step. |
| `GMS_verner_step_avx512.cpp` | AVX-512 Verner Runge–Kutta integration step. |
| `GMS_verner_step_avx512.h` | Declarations for the AVX-512 Verner integration step. |
| `GMS_verner_step_sse.hpp` | Inline/template SSE Verner Runge–Kutta integration step. |
| `misc.h` | Miscellaneous shared kernel declarations and utility definitions. |
