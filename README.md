# RF-EMT — Radio-Frequency Engineering and Modeling Toolkit

<!-- hy-mt2-i18n:start -->

**English** · [中文](./README_zh-CN.md) · [日本語](./README_ja.md) · [Español](./README_es.md)

<!-- hy-mt2-i18n:end -->

[![License: GPL-3.0](https://img.shields.io/badge/License-GPL--3.0-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Top language](https://img.shields.io/github/languages/top/bgin/RF-EMT)](https://github.com/bgin/RF-EMT)
[![Stars](https://img.shields.io/github/stars/bgin/RF-EMT?style=social)](https://github.com/bgin/RF-EMT/stargazers)
[![Forks](https://img.shields.io/github/forks/bgin/RF-EMT?style=social)](https://github.com/bgin/RF-EMT/network/members)

RF-EMT is a research-oriented toolkit for modeling and simulation in radar, telecommunications, and radio-frequency engineering. It brings together electromagnetic models, propagation and scattering utilities, and performance-focused computing techniques to support design-level analysis and integration into larger simulation workflows.

> **Project status:** RF-EMT is under active development. Radar and telecommunications modeling and environmental support libraries are partly implemented; baseband signal processing and optimal receiver modeling are not yet implemented. Feature coverage and build requirements may vary by module.

## Contents

- [Highlights](#highlights)
- [Modeling domains](#modeling-domains)
- [Electromagnetics module](#electromagnetics-module)
- [Performance and toolchains](#performance-and-toolchains)
- [Build and integration](#build-and-integration)
- [Contributing](#contributing)
- [License](#license)

## Highlights

- **RF and radar modeling:** Analytical and simulation components include radar cross-section (RCS) analysis, antenna radiation patterns, propagation, polarization, and system-loss utilities.
- **Electromagnetic and environmental models:** Scattering and dielectric models cover vegetation, hydrometeors, and urban structures, alongside atmospheric and other supporting components.
- **CPU performance paths:** Computational kernels include single- and double-precision implementations using SIMD instruction sets such as SSE, AVX, AVX2, and AVX-512, with OpenMP in selected paths.
- **GPU computing:** The repository includes CUDA kernels and supporting routines for GPU-accelerated computation.
- **Mixed-language scientific computing:** C++, CUDA, and Fortran sources support a range of modeling and numerical requirements.
- **Modular design:** Components can be used as parts of a computational backend or integrated into a larger engineering application or GUI.

## Modeling domains

RF-EMT is organized around these areas:

| Domain | Status |
| --- | --- |
| RF telecommunications | Partly implemented |
| Radar modeling and simulation | Partly implemented |
| Baseband (complex-envelope) processing and optimal receivers | Not yet implemented |
| Supporting libraries, including atmosphere, terrain, and buildings | Partly implemented |

## Electromagnetics module

The [`src/Electromagnetics`](./src/Electromagnetics/) directory is the project's core electromagnetic modeling area. It includes:

- SIMD- and OpenMP-oriented electromagnetic field computation paths
- RCS analysis and scattering models
- Dielectric models for vegetation and other materials
- Urban-structure and tree-scattering components
- Polarization, propagation, and radar-loss utilities
- Fortran T-matrix routines and related scientific code

For an inventory of this module, see [`src/Electromagnetics/README.md`](./src/Electromagnetics/README.md). Atmospheric components and related supporting sources are documented in [`src/Atmosphere/README.md`](./src/Atmosphere/README.md).

## Performance and toolchains

Some kernels are specialized for particular CPU instruction sets or use CUDA. Select compiler options to match the code being built and the hardware on which it will run. Architecture-specific binaries may not run on older or differently configured processors. In particular, `-march=native` produces binaries tuned for the build machine, and fast-math options can change floating-point behavior and numerical results.

The following are examples of commonly used optimization and diagnostic options; they are not a universal RF-EMT build command:

| Toolchain | Example options | Notes |
| --- | --- | --- |
| GCC | `-O3`, `-march=native`, `-ftree-vectorize` | `-ftree-vectorize` is generally enabled by `-O3`. Add ISA-specific flags only when appropriate for the target CPU. |
| Clang/LLVM | `-O3`, `-march=native`, `-mavx2`, `-Rpass=loop-vectorize` | `-mavx2` targets AVX2; use a compatible target and hardware. |
| MSVC | `/O2`, `/arch:AVX2` | Architecture options depend on compiler version and target. |
| Intel oneAPI (ICX) | `-O3`, `-xHost`, `-qopt-report=3` | `-xHost` targets the build host; check compiler documentation for supported options. |

Use relaxed floating-point options such as `-ffast-math` or `/fp:fast` only when their numerical trade-offs are acceptable for your application. CUDA builds also require a compatible NVIDIA toolchain and device. Check the relevant compiler and project configuration for exact requirements.

## Build and integration

RF-EMT is intended as a computational backend for scientific and engineering applications, simulation pipelines, or application frontends. The repository contains modules with different languages and hardware requirements, so build steps and dependencies can differ between components. Consult the relevant project build files and module documentation before compiling; ensure the selected compiler, optional CUDA setup, and CPU instruction-set flags match your target environment.

## Contributing

Contributions are welcome, especially in:

- radar, telecommunications, and radio-frequency engineering
- electromagnetic modeling and numerical methods
- high-performance computing, including SIMD, OpenMP, and CUDA
- scientific software testing, documentation, and build support

Please open an issue to discuss a proposed change, or submit a pull request with a clear description of the problem and solution.

## License

RF-EMT is licensed under the [GNU General Public License, version 3 (GPL-3.0)](https://www.gnu.org/licenses/gpl-3.0).

## Acknowledgments

RF-EMT builds on established engineering and scientific literature, with a focus on realistic and computationally efficient modeling for RF and electromagnetic applications.
