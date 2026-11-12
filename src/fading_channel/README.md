# Fading-channel sources

This directory contains analytical and SIMD-accelerated routines for fading-channel statistics, error probabilities, and related numerical integrals. Files with `.h`/`.cpp` generally separate declarations from implementations; architecture suffixes (`sse`, `avx`, `avx512`) identify optimized variants.

| File | Description |
|---|---|
| `GMS_analytic_bep_sep_ch8.h` | Analytical bit-error and symbol-error probability routines for the chapter 8 fading-channel models. |
| `GMS_analytic_bep_sep_ch8_sse.h` | SSE-optimized analytical bit-error and symbol-error probability routines. |
| `GMS_analytic_sep_bep_ch8_array1d_f32.h` | Single-precision array-based analytical symbol- and bit-error probability calculations. |
| `GMS_analytic_sep_bep_ch8_array1d_f64.h` | Double-precision array-based analytical symbol- and bit-error probability calculations. |
| `GMS_channel_snr_pdf_mgf.h` | Channel SNR probability-density and moment-generating-function utilities. |
| `GMS_compute_functionals_ch4.cpp` | Implementations of fading-channel functionals associated with chapter 4. |
| `GMS_compute_functionals_ch4.h` | Declarations and interfaces for chapter 4 functional calculations. |
| `GMS_compute_functionals_ch5.cpp` | Implementations of fading-channel functionals associated with chapter 5. |
| `GMS_compute_functionals_ch5.h` | Declarations and interfaces for chapter 5 functional calculations. |
| `GMS_fading_spectrum.h` | Fading-channel spectrum definitions and supporting calculations. |
| `GMS_functionals_ch4_sse.cpp` | SSE-optimized implementations of chapter 4 channel functionals. |
| `GMS_functionals_ch4_sse.h` | Declarations for the SSE chapter 4 functional kernels. |
| `GMS_integrands_func_ch4.cpp` | Chapter 4 integrand function implementations used in fading-channel calculations. |
| `GMS_integrands_func_ch4.h` | Declarations for chapter 4 integrand functions. |
| `GMS_integrands_func_ch4_sse.h` | SSE-oriented chapter 4 integrand routines and declarations. |
| `GMS_integrands_func_ch5.cpp` | Chapter 5 integrand function implementations. |
| `GMS_integrands_func_ch5.h` | Declarations for chapter 5 integrand functions. |
| `GMS_integrands_func_ch8.h` | Chapter 8 integrand functions for error-probability calculations. |
| `GMS_integrands_helpers.cpp` | Shared helper implementations used by the integrand routines. |
| `GMS_integrands_helpers.h` | Declarations and utilities shared by fading-channel integrands. |
| `GMS_metaprog_bep_sep_ch8.h` | Generic/template-based chapter 8 bit- and symbol-error probability kernels. |
| `GMS_metaprog_bep_sep_ch8_avx.h` | AVX-specialized template kernels for chapter 8 error probabilities. |
| `GMS_metaprog_bep_sep_ch8_avx512.h` | AVX-512-specialized template kernels for chapter 8 error probabilities. |
| `GMS_metaprog_bep_sep_ch8_sse.h` | SSE-specialized template kernels for chapter 8 error probabilities. |
| `README.md` | This directory guide, listing each file and its role. |
