# HW_profiling folder overview

This directory contains low-level hardware performance profiling utilities for the RF-EMT project. Based on the file names and headers, it focuses on x86 CPU feature detection, performance-monitoring counters (PMCs), timestamp sources (`chrono` and `rdtscp`), SKX/SPR metric preprocessing, and Intel instrumentation helpers.

## File inventory

| File / entry | Type | Short description |
| --- | --- | --- |
| `GMS_avx512_warmup_loops.h` | Header | AVX-512 warm-up code used to stabilize the CPU before timing-sensitive benchmarks. |
| `GMS_chrono_timer.cpp` | C++ source | Implementation of a high-resolution timer built on `std::chrono`. |
| `GMS_chrono_timer.h` | Header | Declaration of the `std::chrono`-based timer class and timing stats helpers. |
| `GMS_cpu_perf_time_series_analysis.h` | Header | Utilities for collecting and analyzing CPU performance time series data. |
| `GMS_cpuid.h` | Header | CPU identification and feature-detection logic for x86/Intel/AMD processors. |
| `GMS_cpuid_x86.c` | C source | Lower-level x86 `CPUID` implementation and processor detection routines. |
| `GMS_fast_pmc_access.h` | Header | Fast access wrappers for reading performance counters and related low-level values. |
| `GMS_hw_perf_macros.h` | Header | Macro-based helpers for hardware counter collection and reporting around code regions. |
| `GMS_perf_collector_chrono.cpp` | C++ source | Implementation of a chrono-based performance collector for timing code sections. |
| `GMS_perf_collector_chrono.h` | Header | Interface for the chrono-based performance collector. |
| `GMS_perf_collector_rdtscp.cpp` | C++ source | Implementation of an `rdtscp`-based performance collector using timestamp counter reads. |
| `GMS_perf_collector_rdtscp.h` | Header | Declarations for the `rdtscp` timing collector. |
| `GMS_pmc_events_descriptive_stats.hpp` | Header | Descriptive statistics helpers for PMU/PMC event data. |
| `GMS_pmc_events_names.h` | Header | Names and definitions for PMC event identifiers. |
| `GMS_postprocess_skx_clk_hw_metrics.hpp` | Header | Post-processing logic for SKX clock and hardware metric data. |
| `GMS_preprocess_skx_clk_hw_metrics.cpp` | C++ source | Preprocessing implementation for Skylake clock-related metrics. |
| `GMS_preprocess_skx_clk_hw_metrics.h` | Header | Header for preprocessing Skylake clock/hardware metrics. |
| `GMS_preprocess_skx_perf_metrics.h` | Header | Header for preprocessing generic Skylake performance metrics. |
| `GMS_rdtscp_timer.cpp` | C++ source | Timing implementation based on the `RDTSCP` instruction. |
| `GMS_rdtscp_timer.h` | Header | Timers using `RDTSCP` and related timing wrappers. |
| `GMS_skx_clk_hw_events_metrics.hpp` | Header | SKX clock-event metric definitions and processing logic. |
| `GMS_skx_cpu_operating_freq_metric.cpp` | C++ source | Measures or calculates CPU operating frequency metrics. |
| `GMS_skx_cpu_operating_freq_metric.h` | Header | Interface for SKX CPU frequency metrics. |
| `GMS_skx_tma_metrics.h` | Header | Intel TMA-related metric definitions for Skylake-class systems. |
| `GMS_spr_metrics_r4.hpp` | Header | SPR (likely Sapphire Rapids) metrics with the `r4` variant. |
| `GMS_spr_metrics_xmm4r4.hpp` | Header | SPR performance metrics for the XMM4R4 vectorized format. |
| `GMS_spr_metrics_ymm8r4.hpp` | Header | SPR performance metrics for the YMM8R4 vectorized format. |
| `GMS_spr_metrics_zmm16r4.hpp` | Header | SPR performance metrics for the ZMM16R4 vectorized format. |
| `GMS_timsac_iface.h` | Header | Interface to TIMSAC-related routines or numerical wrappers. |
| `PerfArch_MSRs.h` | Header | Model-specific register (MSR) definitions for performance architecture monitoring. |
| `SKX_PerfEvtSel_write_masks.h` | Header | Write-mask values used to program SKX performance event select registers. |
| `ittnotify/` | Directory | Intel ITT Notify instrumentation utilities for profiling/tracing support, including a bundled license and static library wrappers. |
| `ittnotify/LICENSE.txt` | License | License text for the bundled `ittnotify` instrumentation package. |
| `ittnotify/disable_warnings.h` | Header | Warning-suppression helpers for integrating Intel ITT Notify cleanly with other code. |
| `ittnotify/ittnotify.h` | Header | Main Intel ITT Notify API header. |
| `ittnotify/ittnotify_config.h` | Header | Configuration definitions for the ITT Notify integration. |
| `ittnotify/ittnotify_static.cpp` | C++ source | Static implementation for ITT Notify functionality. |
| `ittnotify/ittnotify_static.h` | Header | Static library declarations for ITT Notify wrappers. |
| `ittnotify/ittnotify_types.h` | Header | Type definitions used by the ITT Notify API. |
| `ittnotify/legacy/` | Directory | Legacy compatibility code or older instrumentation interfaces for ITT Notify. |

## Overall purpose

This folder is essentially a profiling and hardware-monitoring toolkit for x86 systems. It includes:

- `CPUID`-based processor identification and feature checks
- direct access to hardware counters and MSRs
- timing helpers using both `std::chrono` and `RDTSCP`
- Skylake / Sapphire Rapids performance metric processing
- Intel `ittnotify` instrumentation support for profiling and tracing

This is consistent with a performance-analysis layer intended to measure and characterize low-level CPU behavior during numerical or scientific code execution.
