
# Radio-Frequency Engineering and Modeling Toolkit

<!-- hy-mt2-i18n:start -->

**English** · [中文](./README_zh-CN.md) · [日本語](./README_ja.md) · [Español](./README_es.md)

<!-- hy-mt2-i18n:end -->

This project provides a high-fidelity simulation framework for Radar,Telecommunication and the other Radio-frequency engineering systems. The goal is to achieve realistic, design-level system modeling and simulation, following standards and methodologies found in advanced engineering and technical literature.

## Key Features

- **Comprehensive System Modeling:**  
  Models core components of Radar and Telecommunication systems , including analytical Radar Cross Section (RCS) calculations and antennae radiation patterns.  
- **Performance-Optimized Kernels:**  
  Implements highly optimized algorithms, leveraging Intel Intrinsics (SSE/AVX/AVX2/AVX512) for massive manual vectorization. Compiler-level autovectorization is used for descriptive statistics and profiling.
- **GPGPU Acceleration:**  
  Includes a substantial CUDA codebase (~15,000 lines) covering computational kernels and helper routines for GPU-accelerated simulation.
- **Modular Architecture:**  
  Organized as a collection of standalone modules, each describing distinct modeled components. The library can serve as a computational backend or be integrated with a GUI frontend.
- **Component Scope:**  
  The framework is structured around four main simulation domains:
  1. RF Telecommunication system modeling and simulation (partly implemented).
  2. Radar modeling and simulation (partly implemented).
  3. Baseband (complex envelope) signal processing and optimum receivers modeling (not yet implemented).
  4. Support libraries (e.g. an Atmosphere, terrain and buildings) (partly implemented).
  
## Implementation Overview

- **SIMD Execution Paths:**  
  Hundreds of computational kernels implemented for both double and single precision, focusing on analytical RCS and antenna modeling.
- **CUDA Path:**  
  GPU kernels for large-scale, high-performance computations.
- **Current Status:**  
  The project contains hundreds of SIMD kernels and a comprehensive set of analytical and simulation tools for radar and Telecom system analysis and modeling.

# RF-EMT Fading Channel Directory File List

This markdown file provides an structured overview of the source and header files located within the `src/fading_channel` directory of the **RF-EMT** repository.

| File Name | Extension | Category | Brief Technical Context |
| :--- | :---: | :---: | :--- |
| `GMS_analytic_bep_sep_ch8.h` | `.h` | Header | Analytic Bit Error Probability (BEP) and Symbol Error Probability (SEP) formulations for Chapter 8 context. |
| `GMS_analytic_bep_sep_ch8_sse.h` | `.h` | Header | SSE-optimized implementations of analytic BEP and SEP for Chapter 8. |
| `GMS_analytic_sep_bep_ch8_array1d_f32.h` | `.h` | Header | 1D float32 array configurations for Chapter 8 analytic SEP/BEP models. |
| `GMS_analytic_sep_bep_ch8_array1d_f64.h` | `.h` | Header | 1D float64 (double) array configurations for Chapter 8 analytic SEP/BEP models. |
| `GMS_channel_snr_pdf_mgf.h` | `.h` | Header | Channel Signal-to-Noise Ratio (SNR) Probability Density Function (PDF) and Moment Generating Function (MGF) computations. |
| `GMS_compute_functionals_ch4.h` | `.h` | Header | Header definitions for computing functional math/statistical structures mapped to Chapter 4. |
| `GMS_compute_functionals_ch4.cpp` | `.cpp` | Source | Implementation source for computing mathematical functionals associated with Chapter 4 models. |
| `GMS_compute_functionals_ch5.h` | `.h` | Header | Header definitions for computing functional math/statistical structures mapped to Chapter 5. |
| `GMS_compute_functionals_ch5.cpp` | `.cpp` | Source | Implementation source for computing mathematical functionals associated with Chapter 5 models. |
| `GMS_fading_spectrum.h` | `.h` | Header | Definitions and analytical parameters modeling specialized fading channel spectrum behaviors. |
| `GMS_functionals_ch4_sse.h` | `.h` | Header | Header declaration for SIMD/SSE parallelized implementations of Chapter 4 functionals. |
| `GMS_functionals_ch4_sse.cpp` | `.cpp` | Source | SSE vectorized implementation codebase for Chapter 4 mathematical functionals. |
| `GMS_integrands_func_ch4.h` | `.h` | Header | Integrand function declarations utilizing specific numerical methods for Chapter 4. |
| `GMS_integrands_func_ch4.cpp` | `.cpp` | Source | Implementation of mathematical integrand functions for numerical integration in Chapter 4 models. |
| `GMS_integrands_func_ch4_sse.h` | `.h` | Header | SSE optimized vector parallel versions of the Chapter 4 numerical integrand functions. |
| `GMS_integrands_func_ch5.h` | `.h` | Header | Integrand function declarations utilizing specific numerical methods for Chapter 5. |
| `GMS_integrands_func_ch5.cpp` | `.cpp` | Source | Implementation of mathematical integrand functions for numerical integration in Chapter 5 models. |
| `GMS_integrands_func_ch8.h` | `.h` | Header | Integrand function declarations utilizing specific numerical methods for Chapter 8. |
| `GMS_integrands_helpers.h` | `.h` | Header | Definitions of auxiliary utility routines and wrapper helpers assisting integrand solvers. |
| `GMS_integrands_helpers.cpp` | `.cpp` | Source | Core logic implementations for auxiliary helpers supporting mathematical integration. |


## Usage

This software is intended as a backend computational library for advanced simulation and modeling applications. It can be integrated into larger software environments or connected to graphical user interfaces for visualization and analysis.

## Compiler Infrastructure & Toolchain Configuration Flags

To compile the RF-EMT framework with full hardware optimization enabled, the build system must explicitly target the underlying CPU SIMD instruction sets and enable aggressive vectorization pipelines. Below is the comprehensive guide to configuration flags for the four major C++ compilers:

### 1. GNU Compiler Collection (GCC)
GCC provides fine-grained control over architecture target generation and vectorizer heuristics.
* `-O3`: Activates all high-level optimizations, including aggressive loop vectorization, unrolling, and predictive commoning.
* `-march=native`: Directs the compiler to discover the host CPU topology at compile time and auto-enable all supported SIMD instruction subsets.
* `-mavx512f -mavx512cd -mavx512bw -mavx512dq -mavx512vl`: Explicitly enforces compilation targeting the complete AVX-512 instruction set extension.
* `-ftree-vectorize`: Enables the tree-based auto-vectorization pass (implicitly enabled at `-O3`, but useful for explicit profiling).
* `-ffast-math`: Relaxes strict IEEE 754 compliance to allow algebraic transformations that accelerate floating-point arithmetic (e.g., reciprocal approximations).

### 2. LLVM Clang Compiler
Clang uses an advanced vectorization optimization loop and fully supports native SIMD vector expansions.
* `-O3`: Triggers Clang’s deep pipeline optimizations and multi-pass loop vectorizer.
* `-march=native`: Automatically generates instructions tailored to the compilation host hardware.
* `-mvx2` / `-mavx512f`: Forces code generation for AVX2 or basic AVX-512 foundation blocks.
* `-Rpass=loop-vectorize`: Instructs the optimization pass to output verbose diagnostic remarks regarding which loops were successfully vectorized.
* `-Rpass-missed=loop-vectorize`: Provides detailed feedback on loops that failed vectorization, along with technical reasoning (e.g., pointer aliasing).

### 3. Microsoft Visual C++ (MSVC)
MSVC controls vectorized code generation via specific architecture switches.
* `/O2` / `/Ox`: Enables maximum speed optimizations and full global optimization passes.
* `/arch:AVX2`: Instructs the compiler to emit instructions using AVX2 vector registers.
* `/arch:AVX512`: Directs MSVC to utilize AVX-512 instruction sets for vector operations (available in modern MSVC/Visual Studio toolchains).
* `/Qvec-report:2`: Configures the auto-vectorizer to emit a detailed diagnostic report mapping out exactly which loops were vectorized and which were skipped.
* `/fp:fast`: Relaxes floating-point behavior to unlock substantial speedups during massive matrix manipulation loops.

### 4. Intel oneAPI DPC++/C++ Compiler (ICX)
Intel's modern LLVM-based ICX compiler is uniquely optimized for maximizing SIMD throughput on Intel hardware architectures.
* `-O3`: Instructs the backend optimizer to aggressively vectorize and pipeline loops.
* `-xHost`: Optimizes code execution specifically for the highest SIMD instruction set available on the compilation host machine.
* `-xCORE-AVX512`: Instructs the compiler to output highly efficient AVX-512 instructions designed for Intel Xeon / Core architectures.
* `-qopt-report=3`: Generates a deeply technical, comprehensive optimization report detailing loop transformations, vectorization efficiency, and memory alignments.
* `-fp-model=fast=2`: Maximizes floating-point speed by permitting extensive math modifications and reciprocal approximations.

## Contributing

Contributions are welcome, especially from those with expertise in:
- Radar,Telecommunication and Radioengineering systems
- High-performance computing (SIMD, CUDA)
- Numerical methods and scientific computing

If you are interested, please open an issue or a pull request.

## License

This project is licensed under the terms of the GPLv3 license.

## Acknowledgments

This project builds on a foundation of engineering and technical literature and is driven by a commitment to realistic and efficient system modeling.

