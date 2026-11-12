# Atmosphere sources

This directory contains atmospheric models, weather and cloud microphysics, and supporting numerical utilities. The table covers files directly in this directory; the listed subdirectories contain additional supporting data and model components.

| File or directory | Description |
| --- | --- |
| `CESM_Kernel_License_091118.docx` | License document for CESM-derived kernel code. |
| `GMS_AtmMetricViscosity.h` | Declares the templated atmospheric dynamic-viscosity calculation. |
| `GMS_AtmMetricViscosity.inl` | Implements the temperature-based metric viscosity calculation. |
| `GMS_AtmosphereDefs.h` | Common includes and build definitions for atmospheric C++ code. |
| `GMS_AtmosphereModel76.h` | Declares the U.S. Standard Atmosphere 1976 model interface and data accessors. |
| `GMS_AtmosphereModel76.inl` | Implements the 1976 standard-atmosphere model, including pressure, density, temperature, and sound speed calculations. |
| `GMS_NRLMSISE00_CCOR.h` | Declares NRLMSISE-00 chemistry correction functions. |
| `GMS_NRLMSISE00_CCOR.inl` | Implements the NRLMSISE-00 chemistry correction functions. |
| `GMS_NRLMSISE00_DENSITY.h` | Declares lower-atmosphere temperature and density profile calculations for NRLMSISE-00. |
| `GMS_NRLMSISE00_DENSITY.inl` | Implements NRLMSISE-00 lower-atmosphere temperature and density profiles. |
| `GMS_NRLMSISE00_INOUT.h` | Defines NRLMSISE-00 input, output, switch, and geomagnetic-index structures. |
| `GMS_NRLMSISE00_SCLAHDNET.h` | Declares scale-height and turbopause density-mixing functions for NRLMSISE-00. |
| `GMS_NRLMSISE00_SCLAHDNET.inl` | Implements the NRLMSISE-00 scale-height and density-mixing functions. |
| `GMS_NRLMSISE00_SHARED.h` | Holds shared global state used by the NRLMSISE-00 implementation. |
| `GMS_NRLMSISE00_SPLINE_INTEGRATORS.h` | Declares cubic-spline derivative, interpolation, and integration utilities. |
| `GMS_NRLMSISE00_SPLINE_INTEGRATORS.inl` | Implements the NRLMSISE-00 spline numerical utilities. |
| `GMS_atmos_refraction_state.cpp` | Implements atmospheric refraction state construction and allocation. |
| `GMS_atmos_refraction_state.h` | Defines size parameters and storage state for electromagnetic atmospheric-refraction calculations. |
| `GMS_drag_coeff.cpp` | Implements hydrometeor drag, Reynolds and Davies numbers, and related corrections. |
| `GMS_drag_coeff.h` | Declares hydrometeor drag and flow-parameter calculations. |
| `GMS_hwm14.f90` | Implements the Horizontal Wind Model 2014, including quiet-time and disturbance winds. |
| `GMS_micro_mg2_0.F90` | Implements the Morrison–Gettelman version 2.0 cloud microphysics scheme. |
| `GMS_module_bl_mynn_c_iface.h` | Declares the C interface to the Fortran MYNN boundary-layer module. |
| `GMS_mynn_bl_driver.h` | Provides C++ argument/state wrappers and a call interface for the MYNN boundary-layer driver. |
| `GMS_wind_wrf.f90` | Reads and stores 3-D wind fields from WRF NetCDF output. |
| `GMS_wsm6_driver.f90` | Driver for the WSM6 six-class graupel microphysics scheme. |
| `GMS_wsm6_kernel.F90` | Implements the WSM6 cloud and precipitation microphysics kernel. |
| `micro_mg_cam.F90` | CAM-facing interface that packs inputs, calls MG microphysics, and processes outputs. |
| `micro_mg_utils.F90` | Provides Morrison–Gettelman microphysics constants, particle properties, and process-rate utilities. |
| `module_bl_mynn.F` | Fortran implementation of the MYNN planetary-boundary-layer turbulence and mixing scheme. |
| `shr_kind_mod.F90` | Defines shared Fortran numeric kind constants. |
| `shr_spfn_mod.F90` | Provides shared mathematical special-function utilities, including the gamma function. |
| `wv_sat_methods.F90` | Calculates saturation vapor pressure and humidity over water and ice using selectable methods. |
| `Radiative Transfer/` | Subdirectory for atmospheric radiative-transfer sources. |
| `hwm14-data/` | Subdirectory for coefficient and input data used by HWM14. |
| `micro-mg data/` | Subdirectory for data used by the Morrison–Gettelman microphysics scheme. |

This README summarizes the atmospheric-source modules and supporting data used across the project.
