# Electromagnetics folder overview

This directory contains the electromagnetic-modeling, scattering, and radar-analysis components of the RF-EMT project. The files are heavily vectorized (AVX/AVX2/AVX-512/OpenMP) and include support for electromagnetic field solvers, scattering models, antenna/radar utilities, and Fortran-based transfer-matrix routines.

## Inventory

| File / entry | Type | Short description |
| --- | --- | --- |
| `Antenna/` | Directory | Antenna-related implementation area; not yet expanded in this top-level listing. |
| `EOS/` | Directory | Equation-of-state or electromagnetic/material modeling submodule. |
| `FDTD/` | Directory | Finite-difference time-domain implementations and related routines. |
| `RCS/` | Directory | Radar cross-section (RCS) analysis submodule. |
| `GMS_em_fields_types.h` | Header | Core type definitions used by the electromagnetic field kernels and vectorized code. |
| `GMS_em_fields_xmm2r8.cpp` | C++ source | XMM2R8 implementation of electromagnetic field computations. |
| `GMS_em_fields_xmm2r8.h` | Header | Header for the XMM2R8 electromagnetic-field implementation. |
| `GMS_em_fields_xmm2r8_omp.cpp` | C++ source | OpenMP-accelerated XMM2R8 field solver. |
| `GMS_em_fields_xmm2r8_omp.h` | Header | Header for the OpenMP XMM2R8 implementation. |
| `GMS_em_fields_xmm4r4.cpp` | C++ source | XMM4R4 vectorized electromagnetic field code path. |
| `GMS_em_fields_xmm4r4.h` | Header | Header for the XMM4R4 field implementation. |
| `GMS_em_fields_ymm4r8.cpp` | C++ source | YMM4R8 scalar/vectorized EM field computation implementation. |
| `GMS_em_fields_ymm4r8.h` | Header | Header for the YMM4R8 EM field solver. |
| `GMS_em_fields_ymm4r8_omp.cpp` | C++ source | OpenMP-accelerated YMM4R8 field solver. |
| `GMS_em_fields_ymm4r8_omp.h` | Header | Header for the OpenMP YMM4R8 solver. |
| `GMS_em_fields_ymm8r4.cpp` | C++ source | YMM8R4 implementation of electromagnetic field kernels. |
| `GMS_em_fields_ymm8r4.h` | Header | Header for the YMM8R4 field implementation. |
| `GMS_em_fields_ymm8r4.hpp` | Header | Large inline/templated YMM8R4 EM-field implementation. |
| `GMS_em_fields_ymm8r4_omp.cpp` | C++ source | OpenMP-accelerated YMM8R4 field solver. |
| `GMS_em_fields_ymm8r4_omp.h` | Header | Header for the OpenMP YMM8R4 solver. |
| `GMS_em_fields_ymm8r4_omp.hpp` | Header | OpenMP-optimized inline YMM8R4 implementation. |
| `GMS_em_fields_zmm16r4.cpp` | C++ source | ZMM16R4 vectorized EM field computation implementation. |
| `GMS_em_fields_zmm16r4.h` | Header | Header for the ZMM16R4 field solver. |
| `GMS_em_fields_zmm16r4_omp.cpp` | C++ source | OpenMP-accelerated ZMM16R4 field solver. |
| `GMS_em_fields_zmm16r4_omp.h` | Header | Header for the OpenMP ZMM16R4 solver. |
| `GMS_em_fields_zmm8r8.hpp` | Header | Large zmm8r8 EM field implementation for high-width SIMD execution. |
| `GMS_em_fields_zmm8r8_omp.cpp` | C++ source | OpenMP-backed zmm8r8 solver. |
| `GMS_em_fields_zmm8r8_omp.h` | Header | Header for zmm8r8 OpenMP EM implementation. |
| `GMS_grass_scatterer_AVX512.cpp` | C++ source | AVX-512 grass scatterer simulation code. |
| `GMS_grass_scatterer_AVX512.h` | Header | Data structures and functions for AVX-512 grass scattering. |
| `GMS_grass_scatterers_AVX.cpp` | C++ source | AVX grass-scattering routines for multiple scatterer instances. |
| `GMS_grass_scatterers_AVX.h` | Header | Declarations for AVX grass-scattering logic. |
| `GMS_hydrometeor_ensemble_scattering_AVX.cpp` | C++ source | AVX implementation of hydrometeor ensemble scattering. |
| `GMS_hydrometeor_ensemble_scattering_AVX.h` | Header | Header for AVX hydrometeor scattering routines. |
| `GMS_hydrometeor_ensemble_scattering_AVX512.cpp` | C++ source | AVX-512 hydrometeor ensemble scattering code. |
| `GMS_hydrometeor_ensemble_scattering_AVX512.h` | Header | Header for AVX-512 hydrometeor scattering implementation. |
| `GMS_itm_avx512.cpp` | C++ source | AVX-512 implementation of an ionospheric transmission model (ITM). |
| `GMS_itm_avx512.h` | Header | Header for the AVX-512 ITM implementation. |
| `GMS_matrix_inv_cholesky.cpp` | C++ source | Cholesky-based matrix inversion routines used in EM modeling. |
| `GMS_matrix_inv_cholesky.h` | Header | Declarations for matrix inversion via Cholesky factorization. |
| `GMS_mod_lnf.f90` | Fortran source | Auxiliary Fortran module for log-normal or related linearization routines. |
| `GMS_mod_tmatrix_mps.f90` | Fortran source | Fortran T-matrix / MPS module for scattering calculations. |
| `GMS_naca.cpp` | C++ source | NACA airfoil or profile geometry utilities relevant to scattering/propagation structure. |
| `GMS_naca.hpp` | Header | Header for the NACA-based geometry utilities. |
| `GMS_plant_dielectric_avx512.cpp` | C++ source | AVX-512 plant dielectric property modeling. |
| `GMS_plant_dielectric_avx512.h` | Header | Header for plant dielectric calculations. |
| `GMS_platform_error.hpp` | Header | Platform/compatibility error handling utilities. |
| `GMS_polarization.cpp` | C++ source | Polarization state handling and transformations. |
| `GMS_polarization.h` | Header | Declarations for polarization logic. |
| `GMS_polarization.inl` | Inline template | Inline definitions for polarization operations. |
| `GMS_radar_jamming.hpp` | Header | Radar jamming model and related calculations. |
| `GMS_radar_system_losses.hpp` | Header | System-loss calculations for radar-chain propagation. |
| `GMS_radar_types.h` | Header | Shared radar data structures and enumerations. |
| `GMS_tree_scatterer_AVX.cpp` | C++ source | AVX implementation of tree-scatterer modeling. |
| `GMS_tree_scatterer_AVX.h` | Header | Header for AVX tree-scatterer implementation. |
| `GMS_tree_scatterer_AVX512.cpp` | C++ source | AVX-512 tree-scatterer simulation code. |
| `GMS_tree_scatterer_AVX512.h` | Header | Header for AVX-512 tree-scatterer implementation. |
| `GMS_tree_scatterer_common.cpp` | C++ source | Shared logic for tree-scatterer models and utilities. |
| `GMS_tree_scatterer_common.h` | Header | Common declarations for tree-scatterer routines. |
| `GMS_urban_adt_dyn.h` | Header | Dynamic urban-adjacency/ADT data structures. |
| `GMS_urban_adt_stat.h` | Header | Static urban-adjacency/ADT model definitions. |
| `GMS_urban_adt_stat_avx.h` | Header | AVX-optimized static urban ADT implementation. |
| `GMS_urban_adt_stat_avx512.h` | Header | AVX-512-optimized static urban ADT implementation. |
| `GMS_urban_adt_stat_sse.h` | Header | SSE-optimized static urban ADT implementation. |
| `GMS_urban_model.h` | Header | Core urban electromagnetic model definitions. |
| `GMS_urban_model_avx.h` | Header | AVX-optimized urban model routines. |
| `GMS_urban_model_avx512.h` | Header | AVX-512-optimized urban model routines. |
| `GMS_urban_model_sse.h` | Header | SSE-optimized urban model routines. |
| `lpd.f` | Fortran source | Legacy or supporting Fortran routine for propagation or scattering calculations. |
| `module_tmatrix.f90` | Fortran source | Main Fortran module for T-matrix-based electromagnetic modeling. |
| `tmatrix_mps_np.inc` | Include file | Include support for T-matrix MPS numerical parameters or compile-time configuration. |

## Overall purpose

The `Electromagnetics` folder is the core computational physics section of RF-EMT. It brings together:

- high-performance SIMD electromagnetic field solvers
- scattering and dielectric models for vegetation, hydrometeors, and urban structures
- radar cross-section and radar-loss analysis
- polarization and propagation utilities
- Fortran-based T-matrix and related routines for legacy scientific computations

The naming pattern suggests a strong emphasis on high-performance vectorized physics kernels for propagation and scattering analysis in RF and remote-sensing contexts.
